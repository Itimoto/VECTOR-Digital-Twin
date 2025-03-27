function [doa, doaMat, doaMUSIC] = naive_music(Hest, subcFreq, elemPos, windowSize, thetaRange)
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
    
    % First, define our sizes & establish our Element Position Vector
    [AT, AR, S, K] = size(Hest);

    ulaPos = helper.signedDistances(elemPos); % Get signed distance from each element to the origin
    ulaPos = ulaPos.';                        % Regular transpose for compatibility.

    % Determine the Subcarrier Wavelengths
    c           = physconst('LightSpeed'); 
    subcLambda  = c ./ double(subcFreq);      % Wavelength for each individual subcarrier
    % Compute Array Manifold Vector
    V = zeros(S, AR);
    for s = 1:S
        for ar = 1:AR
            V(s, ar) = -1*2*pi*(ulaPos(ar))/subcLambda(s);
        end
    end

    % Iterate over each snapshot:
    if nargin < 5
        thetaRange = [0 180]; % 0 to 180 degrees (full hemisphere) by default

        if nargin < 4
            windowSize = 1;
        end
    end
    lastValidK = K - windowSize + 1; % When using moving window, this will be the last 'real' snapshot index

    doaMUSIC = zeros(AT, 180*3, S, K);
    doaMat = zeros(AT, AR, S, K);
    for k = 1:lastValidK
        % Now, compute the DOA:
        % For an array centered on the origin, finding the angle off of the
        % array parallel...
        % Remember - Array Manifold Vector is:
        %   V(theta) = exp(-1j*(2pi/lambda)*[y0 y1 ... yN-1]*cos(theta))
        % We model an incoming signal with:
        %   y_at_ar = h_at_ar*X_at = (V(theta))*X_at, so we equate the two.
        %   => cos(theta) = V(theta)/cos(theta) * h_at_ar
        %
        % We are estimating theta for each individual element, so we'll end
        % up with a vector of theta corresponding to each individual element's estimate of DOA
        for at = 1:AT
            for s = 1:S                
                % MUSIC
                % Define possible angles to search through:
                theta = linspace(thetaRange(1), thetaRange(2), 180*3);
    
                % Covariance Matrix:
                %CSI = squeeze(Hest(at, :, s, :)); % Aggregate CSI over both space (AR) and time/snapshots (K)
                % Extract CSI for the current window:
                CSI = squeeze(Hest(at, :, s, k:(k + windowSize - 1)));
                %CSI = squeeze(Hest(at, :, :, k)); % Aggregate CSI over Subcarriers Only? - mixed/noisy results.
                R = (1/size(CSI, 2)) * (CSI * CSI'); % Normalize it

                % Eigen Decomposition
                [eigenVectors, eigenValues] = eig(R);
                [~, idx] = sort(diag(eigenValues), 'descend');
                eigenVectors = eigenVectors(:, idx); % Largest eigenvalue corresponds to Signal, lowest to noise
                noiseSubspace = eigenVectors(:, (1+1):end);

                % MUSIC Spectrum
                musicSpectrum = zeros(length(theta), 1);
                for t = 1:length(theta)
                    V_test = exp(1j * V(s, :) .* cosd(theta(t)).').';
                    musicSpectrum(t) = 1 / abs(V_test' * (noiseSubspace * noiseSubspace') * V_test);
                end

                % Normalize (for Plotting/Testing Purposes)
                musicSpectrum = 10 * log10(musicSpectrum / max(musicSpectrum));
                %musicSpectrum = 10 * log10(musicSpectrum);

                % Extract the maximum
                [~, max_idx] = max(musicSpectrum);
                doaMat(at, :, s, k) = theta(max_idx);
                doaMUSIC(at, :, s, k) = musicSpectrum; % For later analysis
            end
        end
    end

    % Average over doaMat.
    doaAR = doaMat(:, 1, :, 1:lastValidK);  % Duplicates. Just need the first index
    doaAT = mean(doaAR, 1);      % Average between the Transmitting Antennas
    doaS  = mean(doaAT, 3);      % Average along the subcarriers
    doaK  = mean(doaS,  4);      % Average over snapshots

    doa = doaK; % Average it out to yield a number

    % Plot doaMUSIC:
    % Visualize a Single Snapshot
    %{ 
    Implemented in the Plotting Folder.
    for k = 1:1
        at = 1; s = 1;
        figure;
        theta = linspace(0, 180, 180*3);
        plot(theta, doaMUSIC(at, 1:length(theta), s, k));
        xlabel("Theta (deg)");
        ylabel("Likelihood (dB)");
        title("MUSIC Spectrum for each angle Theta (TX: "+at+", SC: "+s+", Snapshot:"+k);
    end
    % Visualize over Snapshots:
    for s = 1:1
        at = 1;
        figure;
        theta = linspace(0, 180, 180*3);
        snapshots = 1:(K - windowSize + 1);
        [X, Y] = meshgrid(snapshots, theta);
        musicSpectrum = squeeze(doaMUSIC(at, 1:length(theta), s, 1:length(snapshots)));
        mesh(X, Y, musicSpectrum);
        xlabel("Snapshot (k)");
        ylabel("DOA (deg)");
        zlabel("Likelihood (dB)");
        title("MUSIC Spectrum over Time for TX: "+at+" and SC: "+s);
    end
    % Visualize over Subcarriers:
    for k = 1:1
        at = 1;
        figure;
        theta = linspace(0, 180, 180*3);
        subcarriers = 1:S;
        [X, Y] = meshgrid(subcarriers, theta);
        musicSpectrum = squeeze(doaMUSIC(at, 1:length(theta), 1:S, k));
        mesh(X, Y, musicSpectrum);
        xlabel("Subcarrier (#)");
        ylabel("DOA (deg)");
        zlabel("Likelihood (dB)");
        title("MUSIC Spectrum over Subcarrier for TX: "+at+" and Snapshot: "+k);
    end
    %}
end