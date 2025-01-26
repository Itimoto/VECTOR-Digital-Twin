function plotDOA_overSnapshots(doaSpectrum, plottitle, thetaRange, windowSize)
    % Plot likelihood pseudospectrum vs. Snapshot (time)
    % Sliders for Transmitting Antenna + Subcarrier
    %
    % DIMENSIONS OF INTEREST
    %   AT  ~ Number of TX Antennas
    %   AR  ~ Number of RX Antennas
    %   S   ~ Number of Subcarriers
    %   K   ~ Number of Snapshots
    %
    %   T   ~ Number of Degrees considered (e.g. linspace(0, 180, T))
    %
    %
    % INPUTS
    % doaSpectrum- CSI Matrix. Dimensions [AT THETA S K]
    % plottitle - String (Plot title)
    % thetaRange- Vector. Contains minimum and maximum values for theta [2]
    % windowSize- (Opt.) The number of snapshots in each averaging window. Plots nicer.
    %

    % First, define our sizes
    if length(size(doaSpectrum) == 4)
        [AT, T, S, K] = size(doaSpectrum);
    else
        [AT, T, S] = size(doaSpectrum);
        K = 1;
    end

    if nargin < 4
        windowSize = 1;
    end

    % Variables for plotting
    theta = linspace(thetaRange(1), thetaRange(2), T);
    snapshots = 1:(K - windowSize + 1);
    [X, Y] = meshgrid(snapshots, theta);

    %%% Set up plot
    fig = figure('Name', plottitle);
    % Set up sliders to select each aspect of the data
    sliderAT = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', AT, 'Value', 1, ...
        'SliderStep', [1/(AT-1), 1/(AT-1)], ...
        'Units', 'normalized', ...
        'Position', [0.1, 0.02, 0.2, 0.05], ...
        'Callback', @(src, ~) updateGraph());
    if AT < 2
        sliderAT.Visible = 'off';
    end

    sliderS = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', S, 'Value', 1, ...
        'SliderStep', [1/(S-1), 1/(S-1)], ...
        'Units', 'normalized', ...
        'Position', [0.7, 0.02, 0.2, 0.05], ...
        'Callback', @(src, ~) updateGraph());
    if S < 2
        sliderS.Visible = 'off';
    end

    % Plot.
    updateGraph();

    % Callback function to update the graph
    function updateGraph()
        % Get slider values
        at = round(sliderAT.Value); % TX Antennas
        s  = round(sliderS.Value);  % Subcarrier

        % Extract + display the selected data
        musicSpectrum = squeeze(doaSpectrum(at, 1:length(theta), s, 1:length(snapshots)));
        mesh(X, Y, musicSpectrum);
        view(0, 90); % View it directly from the top.
        xlabel("Snapshot (k)");
        ylabel("DOA (deg)");
        zlabel("Likelihood (dB)");
        title(plottitle+" | TX: "+at+", SC: "+s);
    end
end