% Asymptotic Maximum Likelihood Estimator
% Based off of Van Trees Vol. 4 (Chapter 8.5: Maximum Likelihood)
function [doa] = naive_maxlikelihood(Hest, elemPos, estRCO, fc, bw)
    % Extracts DOA with the Asymptotic Maximum Likelihood Algorithm
    %
    % Yields DOA relative to Array Parallel (90 = Array Normal, 0 = Array Parallel)
    % (Relative to first element in array) 
    %
    % Assumptions:
    % - ULA (not necessarily uniformly spaced) (spaced in one dimension)
    % - Noise is spatially uncorrelated (no multipath considered)
    % - Condition: S^_f, ml (8.296) may not be non-negative definite
    %
    %
    % DIMENSIONS OF INTEREST
    %   - A: Number of RX Antenna Elements
    %   - K: Number of CSI Packets              ~ assumed to be 1 for the time being
    %   - S: Number of Subcarriers              ~ num_sc
    %
    % Input:
    %   - Hest: CSI used for angle estimation   [A S K]
    %   - elemPos: Element Position             [3 A]
    %   - estRCO: Estimated Radio Chain Offset  [A 1]
    %   - fc: Channel Carrier Frequency (Hz)    (scalar)
    %   - bw: Channel Bandwidth (Hz)            (scalar)
    
    %%% Prepare Data for Analysis
    % Extract the Subcarrier Wavelengths
    c = physconst('LightSpeed');
    S = size(Hest, 2); % Number of Subcarriers corresponds to that position in the Array
    subcarrier_freq = linspace(fc-bw/2, fc+bw/2, S);
    subcarrier_lambda = c ./ subcarrier_freq;

    % Get the Antenna Vector and its length
    A = size(Hest, 1);
    ant_diff = elemPos(:, 2:end) - elemPos(:, 1);     % [3 A-1]
    ant_diff_length = vecnorm(ant_diff);              % [1 A-1] (Distance from each antenna to the next antenna in array)
    ant_diff_length = [0 ant_diff_length];            % [1 A]  (Distance from the first antenna to the first antenna is 0!)

    % Get number of Snapshots, K
    K = size(Hest, 3);

    % Extract CSI Phase
    csi_phase = unwrap(angle(Hest), [], 2); % [A S] Unwrap along the Subcarrier Axis

    % Calculate the Phase Difference
    for j = 2:(size(Hest, 1))
        % Check if the next trace on AP is higher/lower (from same STA antenna)
        if csi_phase(j, 1) > csi_phase(j - 1, 1)
            % Figure out the absolute difference to figure out which mult of 360 should be subtracted
            interDiff = ceil(abs(csi_phase(j, 1) - csi_phase(j - 1, 1)) / (2*pi)); 
            % Subtract 360 to put it all in the same area (to keep deltaPhi constant)
            csi_phase(j, :) = csi_phase(j, :) - (2*pi)*interDiff;
        end
    end

    phase_diff = csi_phase(2:end, :) - csi_phase(1, :) - permute(estRCO(2:end, :), [1 2]);                   % [A S]

    %% Begin Estimation
    doa = zeros(S);   % Output Variable
    for curr_sc = 1:1 % Iterate over all subcarriers. For now, just the one.
        % Consider only the narrow-band case (one subcarrier)
        %Hest = Hest(:, curr_sc, 1); % First snapshot, first subcarrier %<-- This will break down the line. Local variable needed
        Hest = [0; phase_diff(:, curr_sc)]; % This is probably wrong

        % Compute the Covariance Matrix C_X (Average over K snapshots)
        C_X = (1 / K) * (Hest * Hest'); 
        
        % Define the steering matrix for the array
        theta = linspace(-pi/2, pi/2, 180*100); % Az angles from -90 to 90
    
        % Create array manifold matrix (V(Psi)) for all possible angles
        V = zeros(A, length(theta)); % Steering vectors for each theta
        for i = 1:length(theta)
            for ant = 1:A
                % Compute the steering term for each antenna element based
                % on distance
                V(ant, i) = exp(1j * (2*pi*ant_diff_length(ant)) / (subcarrier_lambda(curr_sc) * sin(theta(i))));
                % The steering vector (and subsequent DOA estimation) will
                % be relative to the first element in the array
            end
        end

        % Define the projection matrix P_V (& ortho P_V)
        P_V = zeros(A, A, length(theta));
        P_V_perp = zeros(A, A, length(theta));
        for i = 1:length(theta)
            % Projection matrix onto the range of V(psi)
            P_V(:, :, i) = V(:, i) * (V(:, i)' / (V(:, i)' * V(:, i)));
            P_V_perp = eye(A) - P_V(:, :, i);
        end

        % Now, minimize the cost function to estimate DOA (theta)
        cost = zeros(length(theta), 1);
        for i = 1:length(theta)
            % Compute P_V x C_X x P_V
            P_V_C_X_P_V = P_V(:, :, i) * C_X * P_V(:, :, i)';

            % Compute tr[P_V_perp x C_X] x P_V_perp
            tr_term = trace(P_V_perp * C_X);

            % Compute cost function
            cost(i) = det(P_V_C_X_P_V + (tr_term * P_V_perp)/(1)); % We divide by (N-D) instead of 1
                % Where N is the number of number of params to estimate (1?)
                % And D is the number of incident waves, assumed to be 1?
            cost(i) = -1*log(cost(i));
        end

        % Find the angle (DoA) that minimizes the cost function
        [~, min_idx] = min(cost);

        % Output the estimated DoA
        doa(curr_sc) = theta(min_idx) * 180/pi;
    end

