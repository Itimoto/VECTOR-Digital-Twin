function doa = naive_doa(Hest, spacing)
    % Return Direction of Arrival in Degrees via the ULA Steering Equation
    %   DoA is given relative to the Array Parallel (doa = 90 corresponds to Array Normal)
    %
    % This is an alright approximation, given no noise and no multipath.
    %
    % DIMENSIONS OF INTEREST
    %   SUB ~ Number of Subcarriers
    %   TX  ~ Number of TX Antennas
    %   RX  ~ Number of RX Antennas
    %
    %
    % INPUTS
    % Hest - CSI Matrix. Dimensions [SUB TX RX]
    % spacing - Scalar. Spacing between each Antenna Element

    numSTAant = size(Hest, 2);
    numAPant  = size(Hest, 3);

    phases = rad2deg(unwrap(angle(Hest)));
    pdiff = zeros(size(Hest));

    for i = 1:numSTAant
        for j = 2:(numAPant)
            % Check if the next trace on AP is higher/lower (from same STA antenna)
            if phases(1, i, j) > phases(1, i, j - 1)
                % Figure out the absolute difference to figure out which mult of 360 should be subtracted
                interDiff = ceil(abs(phases(1, i, j) - phases(1, i, j - 1)) / 360); 
                % Subtract 360 to put it all in the same area (to keep deltaPhi constant)
                phases(:, i, j) = phases(:, i, j) - 360*interDiff;
            end
        end
    
        pdiff(:, i, :) = phases(:, i, 1:end) - phases(:, i, [2:end, 1]); % Subtract each element from each following element in ULA
    end

    doa = mean(rad2deg(acos(deg2rad(pdiff)/(2*pi*spacing))));       % Pull Naive DOA from ULA Steering (Array Factor) equation
end