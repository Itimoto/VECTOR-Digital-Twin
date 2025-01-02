%% DOA via Real CSI Data
% Apply data collected by Picoscenes and parsed into a .mat file
% containing:
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
% Done with `naive_doa.m`. Nothing special, just the ULA equation.

% Navigate to the CSI file:
[file, path] = uigetfile('*.mat',...
                    'Select the CSI file [AT AR S K]'); % Stores an `outputMatrix` of the correct dimensions

% Load the MAT file
csiFile = load(fullfile(path, file));
Hest    = csiFile.outputMatrix;
fc      = csiFile.centerFreq;
bw      = csiFile.chanBW;
elemPos = csiFile.elemPos; 

% Pull ULA spacing from user terminal:
if ~exist('elemPos', 'var')
    spacing = input('Enter the spacing between each element in ULA (cm): ');
    spacing = spacing*1e-2; % cm to m
else
    shiftedElemPos = [elemPos(end, :); elemPos(1:end-1, :)]; % Shift everything over
    elemDiff = elemPos - shiftedElemPos; % Get the distance from each element to the next
    meanDiff = mean(elemDiff(2:end, :));
    spacing = vecnorm(meanDiff); % Get the magnitude of the average vector to each next element
end

if ~exist('fc', 'var')
    fc = input('Enter the center frequency for the CSI (MHz): ');
    fc = fc*1e6; 
end

c = physconst('lightspeed');
lambda = c/fc;
spacing = spacing/lambda; % Get spacing in terms of wavelength

% Now, we iterate over the number of snapshots:
doa_per_snapshot = [];%zeros(size(Hest, 4));
for k = 1:size(Hest, 4) % Dimension 4 ~ K (each snapshot)
    currCSI = permute(Hest(:, :, :, k), [3, 1, 2]); % [AT AR S] -> [S AT AR]
    doa_per_snapshot(k, :, :) = doa_lib.naive_doa(currCSI, spacing);
end

disp("Average DOA: " + mean(doa_per_snapshot));

for at = 1:size(Hest, 1)
    figure;
    hold on;
    leg = [];
    % Plot actual DoA estimate
    for ar = 1:(size(Hest, 2) - 1)
        if ar == 2
            %continue;
        end

        plot(abs(doa_per_snapshot(:, at, ar)));
        leg = [leg "TX:"+at+"|DEL:RX:"+ar+"-"+(ar+1)];
    end

    % Now, plot average of DoA estimate
    for ar = 1:(size(Hest, 2) - 1)
        if ar == 2
            %continue;
        end

        % Regular Average
        avg = mean(abs(doa_per_snapshot(:, at, ar)));
        plot(avg*ones(size(doa_per_snapshot, 1)));
        % Moving Average
        %movAvg = movmean((doa_per_snapshot(:, at, ar)), 300); % Window of 10
        %plot(abs(movAvg), '--', 'LineWidth', 1);
        
        leg = [leg "MEAN"+at+"|DEL:RX:"+ar+"-"+(ar+1)];
    end

    title("TX:"+at+" DOA over Time");
    xlabel("Snapshot (K)");
    ylabel("DOA off Array Parallel (Degrees)");    
    legend(leg);
end