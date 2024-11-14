% Based off of Feature Extraction Example from `https://tns.thss.tsinghua.edu.cn/wst/docs/features`
function [doa] = naive_music(Hest, elemPos, estRCO, fc, bw)
    % Extracts DOA with the MUSIC Algorithm
    % Yields DOA relative to Array Parallel (90 = Array Normal, 0 = Array Parallel)
    %
    % DIMENSIONS OF INTEREST
    %   - A: Number of RX Antenna Elements
    %   - T: Number of CSI Packets              ~ assumed to be 1 for the time being
    %   - S: Number of Subcarriers              ~ num_sc
    %
    % Input:
    %   - Hest: CSI used for angle estimation   [A S]
    %   - elemPos: Element Position             [3 A]
    %   - estRCO: Estimated Radio Chain Offset  [A 1]
    %   - fc: Channel Carrier Frequency (Hz)    (scalar)
    %   - bw: Channel Bandwidth (Hz)            (scalar)
    
    %%% Prepare Data for Analysis
    % Extract the Subcarrier Wavelengths
    c = physconst('LightSpeed');
    num_sc = size(Hest, 2); % Number of Subcarriers corresponds to that position in the Array
    subcarrier_freq = linspace(fc-bw/2, fc+bw/2, num_sc);
    subcarrier_lambda = c ./ subcarrier_freq;

    % Extract CSI Phase
    csi_phase = unwrap(angle(Hest), [], 2); % [A S] Unwrap along the Subcarrier Axis

    % Get the Antenna Vector and its length
    ant_diff = elemPos(:, 2:end) - elemPos(:, 1);     % [3 A-1]
    ant_diff_length = vecnorm(ant_diff);              % [1 A-1] (Distance from each antenna to the next antenna in array)
    ant_diff_normalize = ant_diff ./ ant_diff_length; % [3 A-1]

    % Calculate the Phase Difference
    %%%phase_diff = csi_phase(:, 1) - permute(estRCO(2:end, :), [1 2]); % (apply RCO) [A S]
    %{
    phase_diff = csi_phase(2:end, :) - csi_phase(1, :) - permute(estRCO(2:end, :), [1 2]);                   % [A S]
    phase_diff = unwrap(phase_diff, [], 2);
    phase_diff = mod(phase_diff + pi, 2*pi) - pi;
    %}

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
    
    % Broadcasting is performed, get the value of cos(theta) for each packet and each antenna pair
    cos_mat = subcarrier_lambda .* phase_diff ./ (2 .* pi .* permute(ant_diff_length, [2 1]));
    %cos_mat_mean = mean(cos_mat, [2]); % Along the Subcarrier Dimension
    cos_mat_mean = cos_mat; %???

    %%% Perform MUSIC Estimation
    % Construct Signal Covariance MAtrix
    lambda_R = 1e-3;  % Regularization parameter
    R = cos_mat_mean * cos_mat_mean' + lambda_R * eye(size(cos_mat_mean, 1));
    
    % Eigenvalue Decomposition for MUSIC
    [V, D] = eig(R);

    % Sort Eigenvalues in Descending Order
    [~, ind] = sort(diag(D), 'descend');
    V = V(:, ind);

    % Signal Subspace (Corresponding to largest eigenvalues)
    signal_subspace = V(:, 1); % Choosing the number of signal sources (1 source - pick first column)

    % MUSIC Spectrum Calculation
    theta_vals = linspace(-90, 90, 180); % Range of angles from -90 to 90 degrees
    music_spectrum = zeros(size(theta_vals));
    
    for i = 1:length(theta_vals)
        % Calculate the steering vector for each angle (theta)
        % Compute cos(theta) and sin(theta) for each antenna
        steering_vec = exp(1j * 2 * pi * (cosd(theta_vals(i)) * ant_diff_normalize(1, :) + ...
                                         sind(theta_vals(i)) * ant_diff_normalize(2, :)) .* ...
                                         reshape(ant_diff_length, [1, numel(ant_diff_length)]) / c);
        steering_vec = steering_vec(:); % Ensure it's a column vector
        
        % Calculate the MUSIC spectrum (inverse of the projection onto the noise subspace)
        noise_subspace_projection = eye(size(V, 1)) - signal_subspace * signal_subspace';
        
        % Ensure the steering vector and the subspace projection are compatible
        music_spectrum(i) = 1 / abs(steering_vec' * noise_subspace_projection * steering_vec);
    end

    % Find the peak in the MUSIC spectrum
    [~, peak_idx] = max(music_spectrum);
    doa = theta_vals(peak_idx);     % DOA corresponding to the peak

    % Add 90 degrees to convert to relative to Array Parallel
    doa = doa + 90;

    %{
    % Symbolic Nonlinear Optimization is performed:
    syms x y
    % aoa_sol = [x; y; (1-sqrt(x^2 + y^2)];
    aoa_init = [sqrt(1/3); sqrt(1/3); sqrt(1/3)];
    aoa_mat_sol = zeros(3, 1); % packet_num = 1
    options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'none');

    cur_nonlinear_func = @(aoa_sol) ant_diff_normalize' * aoa_sol - cos_mat_mean';
    cur_aoa_sol = lsqnonlin(cur_nonlinear_func, aoa_init, [], [], options);
    aoa_mat_sol = cur_aoa_sol;

    % Yields an [x, y, z] unit vector pointing towards the Signal Source
    doa_mat = aoa_mat_sol ./ vecnorm(aoa_mat_sol); % [3]

    % Now, convert to Theta, Phi coordinates relative to Array Center
    theta = atand(doa_mat(2) ./ doa_mat(1)); % Off Normal
    phi   = acosd(doa_mat(3) ./ sqrt(doa_mat(1).^2 + doa_mat(2).^2 + doa_mat(3).^2)); % Off the Z-Axis

    doa = theta + 90;
    %}
end