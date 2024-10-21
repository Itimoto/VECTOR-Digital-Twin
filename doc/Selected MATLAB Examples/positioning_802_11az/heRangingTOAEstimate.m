function [toa,tau,pow] = heRangingTOAEstimate(chanEstActiveSC,activeFFTIndices,fftLength,sampleRate,numPaths)
%heRangingTOAEstimate Estimates the TOA (Time of Arrival) from a CFR
%
%   [TOA,TAU,POW] = heRangingTOAEstimate(CHANESTACTIVESC,
%   ACTIVEFFTINDICES,FFTLENGTH,SAMPLERATE,NUMPATHS) estimates the TOA
%   from CFR by estimating the CFR's correlation matrix, performing MUSIC
%   super-resolution to estimate its delay profile and detecting the ToA
%   from it.
%
%   TOA is the estimated time-of-arrival in seconds for the packet
%
%   TAU is the estimated time vector for the estimated delay profile
%
%   POW is the estimated delay profile path power vector
%
%   CHANESTACTIVESC is N-by-P array of CFRs with possibly missing
%   subcarriers. N is the active subcarrier CFR FFT length and P is the
%   number of CFR snapshots
%
%   ACTIVEFFTINDICES is a N-by-1 vector of the active FFT indices
%   corresponding to CHANESTACTIVESC
%
%   FFTLENGTH is the full length of the CFR (without any missing
%   subcarriers)
%
%   SAMPLERATE is the sample rate in Hz for the input, CFR
%
%   NUMPATHS is either the estimated or true number of paths in the channel

%   Copyright 2020-2021 The MathWorks, Inc.

narginchk(5,5);
if size(chanEstActiveSC,1)~=numel(activeFFTIndices)
    error('The size of the first dimension of chanEstActiveSC must be equal to activeFFTIndices.');
end

if ~any(fftLength==[256 512 1024 2048])
    error('fftLength must be 256, 512, 1024, or 2048');
end

if ~any(sampleRate==[20 40 80 160]*1e6)
    error('sampleRate must be 20e6, 40e6, 80e6, or 160e6 Hz.');
end

validateattributes(numPaths,{'numeric'},{'scalar','integer','>',0},mfilename,'''numPaths'' value');

subcarrierSpacing = sampleRate/fftLength;
if subcarrierSpacing~=78125
    error('Incorrect combination of fftLength and sampleRate. fftLength must be 256, 512, 1024, or 2048 for sample rate 20e6, 40e6, 80e6, or 160e6.');
end

% Reshape to consider each STS and Rx CFR as an independent
% snapshot
[numActiveSC,numSTS,numRx] = size(chanEstActiveSC);
chanEstActiveSC = reshape(chanEstActiveSC,[numActiveSC,numSTS*numRx]);

% Interpolate across missing subcarriers (MUSIC assumes uniform spacing)
allFFTIndices = (1:fftLength).';
magPart = interp1(activeFFTIndices,abs(chanEstActiveSC),allFFTIndices);
phasePart = interp1(activeFFTIndices,unwrap(angle(chanEstActiveSC)),allFFTIndices);
[realPart,imagPart] = pol2cart(phasePart,magPart);
chanEstWithNan = complex(realPart,imagPart);
chanEst = chanEstWithNan(~isnan(chanEstWithNan(:,1)),:);

% Generate correlation matrix
Rxx = heRangingCorrelationMatrixEstimate(chanEst);

% Estimate multipath delay profile with MUSIC
spatialResolution  = 0.1; % Estimated delay profile spatial resolution (10cm)
reqNumSamples = physconst('LightSpeed')/(spatialResolution*subcarrierSpacing);
reqNumSamples = 2^(round(log2(reqNumSamples))); % Round to closest power of 2 to speed up pmusic
[pow,w] = pmusic(conj(Rxx),numPaths,reqNumSamples,'corr','centered');
tau = w/(2*pi*subcarrierSpacing);% Scale pseudospectrum to time-domain

if ~any(pow==Inf)    
    % Dynamic range for thresholding depending on BW
    % (dynamicRangeArray is determined statistically for different bandwidths,
    % by running longer simulations and estimating the values which return the
    % least ranging error for each bandwidth)
    dynamicRangeArray = [5 6 11 14];
    dynamicRange = dynamicRangeArray(log2(fftLength/256) + 1);
    % Find threshold based on the dynamic range
    threshold = max(pow)/(10.^(dynamicRange/10));

    % Find all path peaks in the delay profile
    [~,idx] = findpeaks(pow,'MinPeakHeight',threshold);
else % pmusic returns Inf in ideal conditions
    idx = find(pow==Inf,1);
end
toa = min(tau(idx)); % First peak is ToA of the DLOS path

end