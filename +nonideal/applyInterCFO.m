function y = applyInterCFO(xArr, permArray, interCFO, fs, t_start)
    % Apply Inter-NIC Center Frequency Offset, to simulate meshed-together
    %  radios.
    % xArr ~ [Time Samples, Antennas]
    % permArray ~ [Num NICs, Antennas-per-NIC] - [[2 3];[1 4]]
    %                   * 2 3 ~ inner two antennas are on the same NIC
    %                   * 1 4 ~ outer two antennas are on another NIC
    % interCFO ~ [Num NICs] - [0, 2500]
    %                   * 0 ~ cfo of 0 on first NIC (2, 3)
    %                   * 2500 ~ cfo of 2500 Hz on second NIC (1, 4)
    % fs ~ (Hz) Sampling Frequency
    % t_start ~ (s) Timestamp indicating when this snapshot begins

    [numNICS, antPerNIC] = size(permArray);
    [numSamples, numAnts] = size(xArr);
    y = zeros(numSamples, numAnts);

    for nic = 1:numNICS
        currCFO = interCFO(nic);

        for ant = 1:antPerNIC
            % For each antenna in the NIC, apply the Frequency Offset:
            antIndex = permArray(nic, ant);
            currX = xArr(:, antIndex);
            y(:, antIndex) = nonideal.freqOffsetSignalWithTime(currX, currCFO, fs, t_start);
        end
    end
end

%% TEST SCRIPT:
%{
fs = 20e6; % Sampling frequency (20MHz)
times = [0, 2.2e-3];
N = (times(2)-times(1))*fs;

xArr = ones(N, 4)*1j*1j*(-1); % Four Antennas
permArray = [[2 3];[1 4]];
interCFO = [-1500, 1000];

figure; hold on;
y_snapshots = []
for i = 1:length(times)
    t_start = times(i);  % Current snapshot timestamp
    
    % Apply frequency shift
    y = nonideal.applyInterCFO(xArr, permArray, interCFO, fs, t_start);
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

disp(["Resulting Carrier Freq. Offset on Ant 1:", median(diff(angle(y_snapshots(:,1))))/(2*pi*(1/fs))]);
disp(["Resulting Carrier Freq. Offset on Ant 2:", median(diff(angle(y_snapshots(:,2))))/(2*pi*(1/fs))]);
disp(["Resulting Carrier Freq. Offset on Ant 3:", median(diff(angle(y_snapshots(:,3))))/(2*pi*(1/fs))]);
disp(["Resulting Carrier Freq. Offset on Ant 4:", median(diff(angle(y_snapshots(:,4))))/(2*pi*(1/fs))]);
%}