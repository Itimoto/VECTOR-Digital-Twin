function y = freqOffsetSignalWithTime(x, f_offset, fs, t_start)
    % freqOffsetSignal Introduces a frequency shift in the signal with an
    %   absolute time reference.
    %
    %   Y = heFreqShiftSignal(X, F_OFFSET, FS) returns the frequency-shifted 
    %   signal Y by mixing X with a complex exponential of frequency F_OFFSET.
    %
    %   X is the input signal (column vector or matrix).
    %   F_OFFSET is the frequency shift in Hz.
    %   FS is the sampling frequency in Hz.
    %   T_START is the timestamp (sec) indicating when this snapshot begins
    %
    %   Y is the frequency-shifted version of X.

    % Time index relative to the snapshot start time
    N = size(x,1);
    n = (0:N-1).';  % Sample indices within the snapshot
    t = t_start + n / fs;  % Absolute time vector
    
    % Compute phase shift with absolute time reference
    phaseShiftOverTime = 2 * pi * f_offset * t;  % Phase accumulation based on absolute time
    
    % Generate complex exponential with tracked phase
    phaseMultiplier = exp(1j * phaseShiftOverTime);
    
    % Apply the frequency shift
    y = x .* phaseMultiplier;
end

%% TEST SCRIPT:
%{
fs = 1e6; % Sampling frequency (1MHz)
cfo = 50e3; % Frequency shift (50kHz)
N = 2048;
times = [0, 2e-3, 4e-3, 6e-3];

x = exp(1j * 2 * pi * 0 * (0:N-1).' / fs); % 0Hz signal - cyclostationary
y_snapshots = [];

figure; hold on;
for i = 1:length(times)
    t_start = times(i);  % Current snapshot timestamp
    
    % Apply frequency shift
    y = nonideal.freqOffsetSignalWithTime(x, cfo, fs, t_start);
    y_snapshots = [y_snapshots, y];

    % Normalize phase to [-pi, pi] range for plotting
    phaseShiftOverTime = mod(angle(y) + pi, 2*pi) - pi;

    % Plot the phase shift for this snapshot
    t = t_start + (0:N-1).' / fs;  % Absolute time vector
    plot(t * 1e3, phaseShiftOverTime * 180/pi, 'DisplayName', ['t_0 = ' num2str(t_start*1e3) ' ms']);
end

grid on;
xlabel('Time (ms)');
ylabel('Phase Shift (degrees)');
title('Phase Shift Over Time Across Snapshots');
legend;

disp(["Resulting Carrier Freq. Offset", median(diff(angle(y_snapshots(:,1))))/(2*pi*(1/fs))]);
%}
