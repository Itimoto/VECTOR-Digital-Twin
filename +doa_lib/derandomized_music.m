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

    if nargin < 5
        windowSize = 3;
    end
    if windowSize < 3 % AT LEAST 3!
        windowSize = 3;
    end
    lastValidK = K - windowSize + 1; % When using moving window, this will be the last 'real' snapshot index

    doaMUSIC = zeros(AT, T, S, K);
    doaMat   = zeros(AT, AR, S, K);
    for k = 1:lastValidK
        % For every frame...

        doaMUSIC_compare = zeros(C, AT, T, S, windowSize);
        for w = 1:(windowSize)
            % We run through a small window into the future...
            
            for combo = 1:C
                % And try every combination for that given frame...
                Hest_combo = Hest(:, COMBINATIONS(combo, :), :, k+w-1); % Rearrange Hest via combo/rearrangement
                Hest_combo = repmat(Hest_combo, [1, 1, 1, 2]);  % Repeat along k to have MUSIC alg work.
    
                [~, ~, doaMUSIC_tmp] = doa_lib.naive_music(Hest_combo, subcFreq, elemPos, 2, thetaRange);
                doaMUSIC_compare(combo, :, :, :, w) = doaMUSIC_tmp(:, :, :, 1); % Put in one of the k's
            end
        end

        % Now, we wish to figure out which sets of frames are the most similar.
        % We choose the starting frame with the least change over each of the frames in the window:
        %{
        similarityMatrix = zeros(C, C, windowSize);
        secondPeakMagnitudes   = zeros(C, windowSize);
        for w = 1:1%(windowSize-1)
            for combo1 = 1:C
                if k == 1
                    % For the first frame, we iterate through each possibility
                    frame1 = reshape(doaMUSIC_compare(combo1, :, :, :, w), [], 1);
                else
                    frame1 = reshape(frame1_prev(1, :, :, :, 1), [], 1);
                end

                for combo2 = 1:C
                    frame2 = reshape(doaMUSIC_compare(combo2, :, :, :, w+1), [], 1);
                    
                    % Determine the similarity between these two matrices:
                    similarityMatrix(combo1, combo2, w) = dot(frame1(:), frame2(:)) / (norm(frame1(:)) * norm(frame2(:))); 
                end

                % Find magnitude of second (normalized) peak, averaged across subcarriers
                avgSecondPeakMag = zeros(1, S);
                for s = 1:S
                    [pkAmp, pkIdx] = findpeaks(squeeze(doaMUSIC_compare(combo1, 1, :, s, w)));
                    [~, idx] = sort(pkAmp, 'descend');
                    if length(idx) > 1
                        avgSecondPeakMag(s) = pkAmp(idx(2));
                    end
                end
                secondPeakMagnitudes(combo1, w) = mean(avgSecondPeakMag);
            end
        end
        %}

        if k == 1
            framePrev = squeeze(doaMUSIC_compare(1, 1, :, :, 1));
            likelyC = 1; % Pick the first one, for now.
        else 
            similarityArr = zeros(1, C);
            for comboCurr = 1:C
                frameCurr = squeeze(doaMUSIC_compare(comboCurr, 1, :, :, 1));

                similarityArr(comboCurr) = mean((framePrev - frameCurr).^2, 'all'); % Get MSE
            end
            [~, likelyC] = min(similarityArr); % Get min MSE along axis (smaller = more similar)
        end
        
        %[~, likelyC] = max(max(max(similarityMatrix, [], 3), [], 2), [], 1);
        % Prep for next iteration:
        %frame1_prev = doaMUSIC_compare(likelyC, :, :, :, :); %reshape(doaMUSIC_compare(likelyC, :, :, :, 1), [], 1);
        %[~, likelyC] = min(min(secondPeakMagnitudes, [], 2), [], 1);

        %for combo = 1:C
        %    doaMUSIC_tmp = squeeze(doaMUSIC_compare(combo, :, :, :, :));
        %    plotting.plotDoA_overSubcarriers(doaMUSIC_tmp, "DOA ESTIMATE - "+num2str(combo), thetaRange);
        %    colorbar;
        %end


        % Assuming non-normalized spectra, we pick the 'most likely' combo:
        % First maximum along AT, then along Subcarriers, then along Theta, and finally between the different combinations.
        %[~, likelyC] = max(max(max(max(doaMUSIC_compare(:, :, :, :, 1), [], 2), [], 4), [], 3), [], 1);
        
        % Deposit the `likelyC` exactly the same way that `naive_music` does.
        theta = linspace(thetaRange(1), thetaRange(2), T);
        for at = 1:AT
            for s = 1:S
                musicSpectrum = doaMUSIC_compare(likelyC, at, :, s, 1);
                [~, max_idx] = max(musicSpectrum);
                doaMat(at, :, s, k) = theta(max_idx);
                doaMUSIC(at, :, s, k) = musicSpectrum;
            end
        end
    end

    % Average over doaMat.
    doaAR = doaMat(:, 1, :, 1:lastValidK);  % Duplicates. Just need the first index
    doaAT = mean(doaAR, 1);      % Average between the Transmitting Antennas
    doaS  = mean(doaAT, 3);      % Average along the subcarriers
    doaK  = mean(doaS,  4);      % Average over snapshots

    doa = doaK; % Average it out to yield a number
end