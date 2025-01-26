function plotCSI(Hest, plottitle, fc, bw)
    % Plot the CSI given
    %
    % DIMENSIONS OF INTEREST
    %   AT  ~ Number of TX Antennas
    %   AR  ~ Number of RX Antennas
    %   S   ~ Number of Subcarriers
    %   K   ~ Number of Snapshots
    %
    %
    % INPUTS
    % Hest      - CSI Matrix. Dimensions [AT AR S] or [AT AR S K]
    % plottitle - String (Plot title)
    % fc        - Scalar. Center/Carrier Frequency (Hz) (used to determine center wavelength)
    % bw        - Scalar. Channel Bandwidth (Hz) (used to determine subcarrier wavelengths)
    %
    % Example usage:
    %   Hest = Hest_Multi(:, :, :, 1); % [AT AR S K] Select the first snapshot
    %   plotting.plotCSI(Hest, "CSI Estimation from Ranging Process. Naive DOA="+num2str(beta_a)+"deg", fc, bw*1e6);
    % (If there are multiple, we add a slider)
    %   plotting.plotCSI(Hest_Multi, "CSI Estimation from Ranging Process. Naive DOA="+num2str(beta_a)+"deg", fc, bw*1e6);

    % First, define our sizes
    if length(size(Hest) == 4)
        [AT, AR, S, K] = size(Hest);
    else
        [AT, AR, S] = size(Hest);
        K = 1;
    end

    if nargin > 2
        % Determine the Subcarrier Frequencies
        c           = physconst('LightSpeed'); 
        bw          = bw;
        subcFreq    = linspace(fc-bw/2, fc+bw/2, S);
    end

    % Define colormap and line styles
    baseColors = lines(AT); % Generate unique base colors for each AT
    %markerPatterns = {'.', '*', 'star', 'square', 'diamond', 'x'}; % Marker patterns for each AT
    lineStyles = {'-', '--', ':', '-.'}; % Line patterns for each AT

    if exist("subcFreq")
        xaxis = subcFreq;
        xtitle = "Subcarrier Index";
    else
        xaxis = linspace(-S/2, S/2, S);
        xtitle = "Subcarrier Frequency (Hz)";
    end

    %%% Set up plot
    fig = figure('Name', plottitle); sgtitle(plottitle);
    rows = 2;
    % Axes for plots
    axMag = subplot(rows, 1, 1, 'Parent', fig);
    axPhase = subplot(rows, 1, 2, 'Parent', fig);

    % Slider UI to select Snapshot K
    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', K, 'Value', 1, ...
        'SliderStep', [1/(K-1), 1/(K-1)], ...
        'Units', 'normalized', ...
        'Position', [0.3, 0.02, 0.4, 0.05], ...
        'Callback', @(src, event) updatePlot(src.Value));
    
    if K < 2 % Disable slider if we don't need it.
        slider.Visible = 'off'; 
    end

    % Initial plot
    updatePlot(slider.Value);

    % Nested function to update plots
    function updatePlot(kSlider)
        % Clear previous plots
        cla(axMag);
        cla(axPhase);

        % Parse input
        if kSlider < 1
            k = 1;
        elseif kSlider > K
            k = K;
        else
            k = floor(kSlider);
        end
        
        % Plot Magnitude
        subplot(rows, 1, 1);
        hold(axMag, 'on'); leg = [];
        for at = 1:AT
            for ar = 1:AR
                color = baseColors(at, :) * (0.4 + 0.7 * (ar - 1) / (AR - 1)); % Adjust brightness
                plot(axMag, xaxis, mag2db(abs(squeeze(Hest(at, ar, :, k)))), ...
                    'Color', color, 'LineStyle', lineStyles{mod(at-1, length(lineStyles)) + 1}, 'LineWidth', 1);
                    %'Color', color, 'Marker', markerPatterns{mod(at-1, length(markerPatterns)) + 1}, 'LineWidth', 1);
                leg = [leg, "AT-"+at+"|AR-"+ar];
            end
        end
        title(axMag, "CSI Magnitude"); xlabel(axMag, xtitle); ylabel(axMag, "Magnitude (dB)"); grid(axMag, 'on'); legend(axMag, leg);
    
        % Plot Phase
        hold(axPhase, 'on'); leg = [];
        for at = 1:AT
            for ar = 1:AR
                color = baseColors(at, :) * (0.4 + 0.7 * (ar - 1) / (AR - 1)); % Adjust brightness
                % The first term (0.4) sets minimum brightness, second term (0.6) sets range of brightness
                plot(axPhase, xaxis, rad2deg(angle(squeeze(Hest(at, ar, :, k)))), ...
                    'Color', color, 'LineStyle', lineStyles{mod(at-1, length(lineStyles)) + 1}, 'LineWidth', 1);
                    %'Color', color, 'Marker', markerPatterns{mod(at-1, length(markerPatterns)) + 1}, 'LineWidth', 1);
                leg = [leg, "AT-"+at+"|AR-"+ar];
            end
        end
        title(axPhase, "CSI Phase"); xlabel(axPhase, xtitle); ylabel(axPhase, "Phase (deg)"); grid(axPhase, 'on'); legend(axPhase, leg);
    end
end