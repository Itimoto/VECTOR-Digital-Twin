function [doa, doaMat, doaMUSIC] = naive_doa_plus(Hest, fc, bw, elemPos)
    % Return Direction of Arrival in Degrees via the ULA Steering Equation
    %   DoA is given relative to the Array Parallel (doa = 90 corresponds to Array Normal)
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
    % fc        - Scalar. Center/Carrier Frequency (used to determine center wavelength)
    % bw        - Scalar. Channel Bandwidth (used to determine subcarrier wavelengths)
    % elemPos   - Positions for each RX antenna [AR 3] (Assume a ULA)
    %
    
    % First, define our sizes & establish our Element Position Vector
    [AT, AR, S, K] = size(Hest);

    ulaPos = helper.signedDistances(elemPos); % Get signed distance from each element to the origin
    ulaPos = ulaPos.';                        % Regular transpose for compatibility.

    % Determine the Subcarrier Wavelengths
    c           = physconst('LightSpeed'); 
    subcFreq    = linspace(fc-bw/2, fc+bw/2, S);
    subcLambda  = c ./ subcFreq;                % Wavelength for each individual subcarrier
    % Compute Array Manifold Vector
    V = zeros(S, AR);
    for s = 1:S
        for ar = 1:AR
            V(s, ar) = -1*2*pi*(ulaPos(ar))/subcLambda(s);
        end
    end

    % Iterate over each snapshot:
    windowSize = 200; % Number of snapshots in each window
    doaMUSIC = zeros(AT, 180*3, S, K);
    doaMat = zeros(AT, AR, S, K);
    for k = 1:(K - windowSize + 1)
    %for k = 1:K
    %for k = 1:1
        % First, unwrap phases and place all under the first element in ULA
        phases = rad2deg(unwrap(angle(Hest(:, :, :, k))));

        for at = 1:AT
            for ar = 2:AR
                % Check if the next trace on AP is higher/lower (from same STA antenna)
                if phases(at, ar, 1) > phases(at, ar - 1, 1)
                    % Figure out the absolute difference to figure out which mult of 360 should be subtracted
                    interDiff = ceil(abs(phases(at, ar, 1) - phases(at, ar - 1, 1)) / 360);
                    % Subtract 360 to put it in the same area (keep delta Phi constant)
                    phases(at, ar, :) = phases(at, ar, :) - 360*interDiff;
                end
            end
        end

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
                % Steering Vector Method (Deterministic) (Does not work)
                % Array Manifold Vector, 'unique' to each Subcarrier (due to electrical stretching/shrinking by wavelength)
                % Note - We're only working with angles -- so, we eliminate the exponential/complex term
                % now computed once V = -1*(2*pi)*(ulaPos/subcLambda(s)); % Everything but the cos(theta) term
                %V = -2*pi*(ulaPos) ./ subcLambda(s);
                %V = V - V(1);   % Shift Phase Reference to First Element              
                %localH = rad2deg(unwrap(angle(Hest(at, :, s, k))));
                %localH = phases(at, :, s);   % h11, h12, h13, h14 for current at = 1 and AR = 4 RX antennas
                %if mod(size(localH, 2), 2) % Odd N
                %    phaseAtOrigin = localH((size(localH, 2) + 1) / 2); 
                %else % Even N
                %    lowIndex = floor(size(localH, 2)/2);
                %    phaseAtOrigin = mean(localH); % Average them out to get phase at origin
                %end
                %doaMat(at, :, s, k) = rad2deg(acos(deg2rad(localH - phaseAtOrigin)./V)); % Steering Equation
                % ^^^ That didn't work for some reason

                %{
                % Correlation Method - Only on a Per-Element Basis 
                % (No Spatial Correlation -- Doesn't Take Advantage of Multiple Antennas)
                % Define possible angles to search through:
                theta = linspace(0, 180, 180*3);

                % Compute every steering matrix possible
                currVec = exp(1j * V(s, :) .* cosd(theta).');
                
                % Compute the correlation matrix (inner product) for all theta values at once
                correlations = abs(currVec * Hest(at, :, s, k)'); 

                [~, max_idx] = max(correlations);
                doaMat(at, :, s, k) = theta(max_idx);
                %}
                % MUSIC
                % Define possible angles to search through:
                theta = linspace(0, 180, 180*3);

                % Covariance Matrix:
                %CSI = squeeze(Hest(at, :, s, :)); % Aggregate CSI over both space (AR) and time/snapshots (K)
                % Extract CSI for the current window:
                CSI = squeeze(Hest(at, :, s, k:(k + windowSize - 1)));
                R = CSI * CSI';

                % Eigen Decomposition
                [eigenVectors, eigenValues] = eig(R);
                [~, idx] = sort(diag(eigenValues), 'descend');
                eigenVectors = eigenVectors(:, idx); % Largest eigenvalue corresponds to Signal, lowest to noise
                noiseSubspace = eigenVectors(:, 2:end);

                % MUSIC Spectrum
                musicSpectrum = zeros(length(theta), 1);
                for t = 1:length(theta)
                    V_test = exp(1j * V(s, :) .* cosd(theta(t)).').';
                    musicSpectrum(t) = 1 / abs(V_test' * (noiseSubspace * noiseSubspace') * V_test);
                end

                % Normalize (for Plotting/Testing Purposes)
                musicSpectrum = 10 * log10(musicSpectrum / max(musicSpectrum));

                % Extract the maximum
                [~, max_idx] = max(musicSpectrum);
                doaMat(at, :, s, k) = theta(max_idx);
                doaMUSIC(at, :, s, k) = musicSpectrum; % For later analysis

                %{
                % Another Correlation Method - Based on Matlab's fminunc
                % (Also no Spatial Correlation?)
                theta0 = 90;
                objectiveFcn = @(theta) -abs(exp(1j * V .* cosd(theta).' * Hest(at,:, s, k)'));
                
                options = optimset('Display', 'off'); % Suppress output
                thetaOpt = fminunc(objectiveFcn, theta0, options);
                doaMat(at, :, s, k) = thetaOpt;
                %}
            end
        end
    end

    doaAR = mean(doaMat, 2);    % Average along the AR (Receive Antennas) dimensions
    doaS  = mean(doaAR, 3);     % Average along the subcarriers
    doaAT = mean(doaS, 1);      % Average between transmitters
    doaK  = mean(doaAT, 4);     % Average through time

    doa = doaK; % Average it out to yield a number

    % Plot doaMUSIC:
    for s = 1:1
        at = 1; %s = 26;
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
end