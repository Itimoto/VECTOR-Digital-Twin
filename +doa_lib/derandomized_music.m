function [doa, doaMat, doaMUSIC] = derandomized_music(Hest, subcFreq, elemPos, windowSize, thetaRange)
    % Wrapper for `naive_music`. Assumes Hest has 'randomized' data, with
    %   Receive Antennas Flipped randomly.
    %
    % We assume the setup for the UAkron Testbed. Namely:
    %   - 2 NICs in use (21 and 23). Each NIC has a (M)ain and an Au(X) antenna. When
    %       collecting CSI, each NIC randomly flips the M and X antennas.
    %   - The NICs are arranged like so:
    %       (23X) (21X) (21M) (23M)
    %   - X and M can switch randomly on each. There are 2*2 reconfigurations:
    %       [0 1 2 3] [0 2 1 3] [3 2 1 0] [3 1 2 0]
    % We will try to use every possible combination, then apply it to each
    %   frame to find the pseudospectra. We will compare the pseudospectra
    %   across frames to estimate the 'right' configuration for the given
    %   frame.
    %
    % Return Direction of Arrival in Degrees via the ULA Steering Equation
    %   DoA is given relative to the Array Parallel (doa = 90 corresponds to Array Normal)
    % NOTE - THIS ASSUMES ONLY ONE SIGNAL OF INTEREST! Theoretically can be
    %   adjusted for >1 signal source
    %
    % This is an alright approximation, given no noise and no multipath.
    %
    % DIMENSIONS OF INTEREST
    %   AT  ~ Number of TX Antennas
    %   AR  ~ Number of RX Antennas
    %   S   ~ Number of Subcarriers
    %   K   ~ Number of Snapshots
    %
    %
    % INPUTS
    % Hest      - CSI Matrix. Dimensions [AT AR S K]
    % subcFreq  - Array. Dimensions [1 S]. Subcarrier Frequencies for each S (Hz)
    % elemPos   - Positions for each RX antenna [AR 3] (Assume a ULA)
    % windowSize- The number of snapshots in each averaging window
    % thetaRange- [minAngle maxAngle] (deg) - Angles off of Array Parallel to run MUSIC on
    %
    
    COMBINATIONS = [[1 2 3 4], 
                    [1 3 2 4],
                    [4 3 2 1],
                    [4 2 3 1]];
    [~, C] = size(COMBINATIONS); % C corresponds to the number of combinations we're running

    % First, define our sizes & establish our Element Position Vector
    [AT, AR, S, K] = size(Hest);
    T = 180*3; % Number of angles (T)heta sampled

    doaMUSIC = zeros(AT, T, S, K);
    doaMat   = zeros(AT, 1, S, K);
    for k = 1:(K - 1) % One away from the end
        % For every frame...

        doaMUSIC_compare = zeros(C, AT, T, S, 2);
        for w = 1:2
            % We operate one frame into the future.
            
            for combo = 1:C
                % And try every combination for that given frame...
                Hest_combo = Hest(:, COMBINATIONS(combo, :), :, k+w-1); % Rearrange Hest via combo/rearrangement
                Hest_combo = repmat(Hest_combo, [1, 1, 1, 2]);  % Repeat along k to have MUSIC alg work.
    
                [~, ~, doaMUSIC_tmp] = doa_lib.naive_music(Hest_combo, subcFreq, elemPos, 2, thetaRange);
                doaMUSIC_compare(combo, :, :, :, w) = doaMUSIC_tmp(:, :, :, 1); % Put in one of the k's
            end
        end

        if k == 1
            %framePrev = squeeze(doaMUSIC_compare(1, 1, :, :, 1));
            framePrev = permute(doaMUSIC_compare(1, 1, :, :, 1), [3 4 1 2 5]); % Remove 1, 2, & 5 dimensions
            likelyC = 1; % Pick the first one, for now.
        else 
            similarityArr = zeros(1, C);
            for comboCurr = 1:C
                %frameCurr = squeeze(doaMUSIC_compare(comboCurr, 1, :, :, 1));
                frameCurr = permute(doaMUSIC_compare(comboCurr, 1, :, :, 1), [3 4 1 2 5]); % Remove 1, 2, & 5 dimensions
            
                similarityArr(comboCurr) = mean((framePrev - frameCurr).^2, 'all'); % Get MSE
            end
            [~, likelyC] = min(similarityArr); % Get min MSE along axis (smaller = more similar)
        end
        
        % Deposit the `likelyC` exactly the same way that `naive_music` does.
        theta = linspace(thetaRange(1), thetaRange(2), T);
        %musicSpectrum = reshape(doaMUSIC_compare(likelyC, :, :, :, 1), [], S);
        %musicSpectrum = squeeze(doaMUSIC_compare(likelyC, :, :, :, 1));
        musicSpectrum = permute(doaMUSIC_compare(likelyC, :, :, :, 1), [2 3 4 1 5]); % Effectively truncate 1st & 5th dim.
        [~, max_idx] = max(musicSpectrum, [], 2);
        doaMat(:, 1, :, k) = repmat(theta(max_idx), [1, size(doaMat, 2), 1]);
        doaMUSIC(:, :, :, k) = musicSpectrum;
    end

    % Average over doaMat.
    doaAR = doaMat(:, 1, :, 1:(K-1));  % Duplicates. Just need the first index
    doaAT = mean(doaAR, 1);      % Average between the Transmitting Antennas
    doaS  = mean(doaAT, 3);      % Average along the subcarriers
    doaK  = mean(doaS,  4);      % Average over snapshots

    doa = doaK; % Average it out to yield a number
end