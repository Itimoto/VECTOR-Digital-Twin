function [txSig, srcBits] = helperIndoorRayTracingWaveformGen(numCodewordsPerFrame, cfgLDPCEnc, modOrder, ofdmMod)
%HELPERINDOORRAYTRACINGWAVEFORMGEN Generate waveform

%   Copyright 2020-2021 The MathWorks, Inc. 

% LDPC encoding
msgLen = cfgLDPCEnc.NumInformationBits;
codewordLen = cfgLDPCEnc.BlockLength;
srcBits = randi([0, 1], numCodewordsPerFrame*msgLen, 1);
encBits = zeros(numCodewordsPerFrame*codewordLen, 1);
for i = 1:numCodewordsPerFrame
    encBits((i-1)*codewordLen+1:i*codewordLen,:) = ...
        ldpcEncode(srcBits((i-1)*msgLen+1:i*msgLen,:), cfgLDPCEnc);
end

% QAM modulation
encSym = qammod(encBits, modOrder, ...
    'InputType', 'bit', 'UnitAveragePower', true);

% OFDM modulation using dummy QPSK pilot signals. 
ofdmInfo = info(ofdmMod);
ofdmData = reshape(encSym, ofdmInfo.DataInputSize);
ofdmPilot = qammod( ...
    randi([0, 3], ofdmInfo.PilotInputSize), 4, 'UnitAveragePower', true);
txSig = ofdmMod(ofdmData, ofdmPilot);

% Scaling to normalize the signal power
scalingFactor = ofdmMod.FFTLength/ ... 
    sqrt(ofdmMod.FFTLength-sum(ofdmMod.NumGuardBandCarriers)-1); 
txSig = scalingFactor * txSig; 

% [EOF]