function y = swapTracesRandom(xArr, permArray, minSwitchTime, fs, t_start)
    % Sometimes, NIC antennas swap themselves internally. We simulate that
    %   behavior here.
    % Swap the traces in xArr at random.
    % Every `minSwitchTime`, flip a weighted coin to swap each pair of
    % antennas (as specified in `permArray`)
    if minSwitchTime == inf
        y = xArr;
        return; % Skip if necessary
    end


    [numNICS, antPerNIC] = size(permArray);
    [numSamples, numAnts] = size(xArr);
    
    % Calculate total time in seconds.
    totalTime = t_start + numSamples/fs;
    % Determine number of complete segments in the trace.
    numFlips = floor(totalTime/minSwitchTime);
    % Compute number of samples corresponding to one segment.
    segmentLen = round(minSwitchTime * fs);
    
    % Initialize output with the input data.
    y = xArr;
    
    % Process each segment.
    for kp = 1:numFlips
        % Define the start and end indices for the current segment.
        startIdx = (kp-1)*segmentLen + 1;
        endIdx = min(kp*segmentLen, numSamples);
        
        % For each NIC, decide whether to swap based on a coin flip.
        for nic = 1:numNICS
            if rand >= 0.5  % 50% chance to swap for this NIC in the segment.
                % Get the antenna indices to swap for this NIC.
                antIndex1 = permArray(nic, 1);
                antIndex2 = permArray(nic, 2);
                % Swap the segments for the two antennas.
                temp = y(startIdx:endIdx, antIndex1);
                y(startIdx:endIdx, antIndex1) = y(startIdx:endIdx, antIndex2);
                y(startIdx:endIdx, antIndex2) = temp;
            end
        end
    end
end

%% TEST SCRIPT:
%{
fs = 20e6;
times = [0, 2.2e-3];
N = (times(2) - times(1))*fs;

xArr = ones(N,4);
permArray = [[2 3];[1 4]];
minSwitchTime = 0.1e-3;
xArr(:, 2) = xArr(:, 1) * exp(1j*pi/4);
xArr(:, 3) = xArr(:, 1) * exp(1j*pi/2);
xArr(:, 4) = xArr(:, 1) * exp(1j*3*pi/4);

figure; hold on;
y_snapshots = [];
for i = 1:length(times)
    t_start = times(i);  % Current snapshot timestamp
    
    % Apply random swap
    y = nonideal.swapTracesRandom(xArr, permArray, minSwitchTime, fs, t_start);
    y_snapshots = [y_snapshots, y];

    % Plot the phase shift for this snapshot
    t = t_start + (0:N-1).' / fs;  % Absolute time vector
    plot(t * 1e3, angle(y(:,1)) * 180/pi, 'DisplayName', ['t_0 = ' num2str(t_start*1e3) ' ms - 1']);
    plot(t * 1e3, angle(y(:,2)) * 180/pi, 'DisplayName', ['t_0 = ' num2str(t_start*1e3) ' ms - 2']);
    plot(t * 1e3, angle(y(:,3)) * 180/pi, 'DisplayName', ['t_0 = ' num2str(t_start*1e3) ' ms - 3']);
    plot(t * 1e3, angle(y(:,4)) * 180/pi, 'DisplayName', ['t_0 = ' num2str(t_start*1e3) ' ms - 4']);
end
grid on;
xlabel('Time (ms)');
ylabel('Phase Shift (degrees)');
title('Phase Shift Over Time Across Snapshots');
legend;
%}