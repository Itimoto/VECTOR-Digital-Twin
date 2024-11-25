function pos = signedDistances(elemPos)
    % Computes signed distances of points from the center of line of best
    %   fit in 3D space
    %
    % INPUT:
    %   elemPos - Positions for each element size [N 3]
    %
    % OUTPUT:
    %   pos     - [N 1] vector of signed distancs from each point to the
    %               center along the line of best fit
    %
    % Example:
    %
    % elemPos = [
    %       0 -0.0525 0;
    %       0 -0.0175 0;
    %       0  0.0175 0;
    %       0  0.0525 0
    %   ]; % Returns:
    % pos = [
    %       -0.0525
    %       -0.0175
    %        0.0175
    %        0.0525
    %   ];
    
    % Compute mean position (center of the line)
    center = mean(elemPos, 1);

    % Center positions by subtracting the mean
    centeredPos = elemPos - center;

    % Find Direction of Best-Fit Line
    [~, ~, V] = svd(centeredPos, 'econ');

    % First column of V is the direction of the line of best fit,
    %   project each point onto the line of best fit
    pos = centeredPos * V(:, 1);
end