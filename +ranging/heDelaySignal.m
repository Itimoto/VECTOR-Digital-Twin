function y = heDelaySignal(x,delay)
%heDelaySignal Introduces time delay in the signal
%
%   Y = heDelaySignal(X,DELAY) returns the delayed signal Y by applying
%   DELAY to the input signal X. 
%
%   X represents the input signal that needs to be delayed, which can be a
%   column vector or a matrix.
%
%   DELAY is an integer or non-integer scalar value (in samples).
%
%   Y is a delayed version of X which can be a column vector or a matrix.
%   When X is a matrix with multiple columns, each column of X is delayed
%   by DELAY and the resulting signal will be stored in corresponding
%   column of Y.

%   Copyright 2020 The MathWorks, Inc.

% Initialization
inputLength = size(x,1);     % Input signal length
delayInt = round(delay);     % Integer delay
delayFrac = delay-delayInt;  % Fractional delay
outputLength = inputLength+delayInt; % Output signal length
y = complex(zeros(outputLength,size(x,2)));

% Perform delay operation
if delayFrac
    nfft = 2^nextpow2(outputLength);
    binStart = floor(nfft/2);
    % Notice the FFT bins must belong to [-pi, pi].
    fftBin = 2*pi*ifftshift(((0:nfft-1)-binStart).')/nfft;
    tmpxd = fft(x,nfft);
    tmpxd = ifft(tmpxd.*exp(-1i*delay*fftBin));

    newStart = delayInt+1;
    y(newStart:outputLength,:) = tmpxd(newStart:outputLength,:);
else
    % Integer sample shift
    newStart = delayInt+1;
    y(newStart:outputLength,:) = x;
end
end
