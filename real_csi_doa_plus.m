%% DOA via Real CSI Data
% Apply data collected by Picoscenes and parsed into a .mat file containing:
%
% outputMatrix: Matrix of dimensions [AT AR S K] where:
%   AT: Number of Transmit Antennas
%   AR: Number of Receive Antennas
%   S : Number of Subcarriers
%   K : Number of CSI Snapshots
%
% centerFreq: Center/Carrier Frequency of CSI
% chanBW: Channel Bandwidth
% elemPos: Element Position Matrix of dimensions [3 AR] where:
%   AR: Number of Receive Antennas
%   3 : [X Y Z] position of each Receive Element corresponding to each AR position
%
% Done with `naive_doa_plus.m`. Nothing special, just the ULA equation.
%  This has the added benefit of considering the Array Geometry by Subcarrier

% Navigate to the CSI file:
[file, path] = uigetfile('*.mat',...
                    'Select the CSI file [AT AR S K]'); % Stores an `outputMatrix` of the correct dimensions

% Load the MAT file
csiFile = load(fullfile(path, file));
Hest    = csiFile.outputMatrix;
fc      = csiFile.centerFreq;
bw      = csiFile.chanBW;
elemPos = csiFile.elemPos; 

% Now, perform DOA
[doa, doaMat] = doa_lib.naive_doa_plus(Hest, fc, bw, elemPos);
doaK = mean(doaMat, 3);

disp("Average DOA: " + doa);

for at = 1:size(Hest, 1)
    figure;
    hold on;
    leg = [];
    % Plot actual DoA estimate
    for ar = 1:(size(Hest, 2))
        if ar == 2
            %continue;
        end

        toPlot = abs(squeeze(doaK(at, ar, 1, :)));
        plot(toPlot);
        leg = [leg "TX:"+at+"|DEL:RX:"+ar+"-"+(ar+1)];
    end

    % Now, plot average of DoA estimate
    %{
    for ar = 1:(size(Hest, 2))
        if ar == 2
            %continue;
        end

        % Regular Temporal Average
        toPlot = abs(squeeze(doaK(at, ar, 1, :)));
        avg = mean(toPlot);
        plot(avg*ones(size(doaMat, 1)), '-', 'LineWidth', 1);
        % Moving Average
        %movAvg = movmean(toPlot, 300); % Window of 10
        %plot(abs(movAvg), '--', 'LineWidth', 1);
        
        leg = [leg "MEAN"+at+"|DEL:RX:"+ar+"-"+(ar+1)];
    end
    %}

    title("TX:"+at+" DOA over Time");
    xlabel("Snapshot (K)");
    ylabel("DOA off Array Parallel (Degrees)");    
    legend(leg);
end