function [decBits, eqSym] = helperIndoorRayTracingRxProcessing(rxSig, CIR, chanInfo, ...
    cfgLDPCDec, modOrder, ofdmDemod, snr)
%HELPERINDOORRAYTRACINGRXPROCESSING Perform receiver processing

%   Copyright 2020-2022 The MathWorks, Inc. 

% Retrieve OFDM parameters
fftLen = ofdmDemod.FFTLength;
cpLen = ofdmDemod.CyclicPrefixLength;
numGuardBandCarriers = ofdmDemod.NumGuardBandCarriers;
pilotCarrierIdx = ofdmDemod.PilotCarrierIndices;
numOFDMSymbols = ofdmDemod.NumSymbols;
dataCarrierIdx = setdiff( ...
    numGuardBandCarriers(1)+1:fftLen-numGuardBandCarriers(2), ...
    [pilotCarrierIdx; fftLen/2+1]);

% Perfect channel estimation
chanDelay = channelDelay(CIR, chanInfo.ChannelFilterCoefficients);
chanEst = ofdmChannelResponse(CIR,chanInfo.ChannelFilterCoefficients, ...
    fftLen, cpLen, dataCarrierIdx, chanDelay);

% Revert the power normalization at Tx
scalingFactor = sqrt(fftLen-sum(numGuardBandCarriers)-1)/fftLen;
rxSig = scalingFactor * rxSig; 

% OFDM demodulation
ofdmInfo = info(ofdmDemod);
rxSig = rxSig(chanDelay+(1:ofdmInfo.InputSize(1)),:);
rxOFDM = ofdmDemod(rxSig);

% Calculate noise variance
numTx = size(chanEst, 3);
numRx = size(chanEst, 4);
sigPower = numTx/numRx;
nVar = (scalingFactor^2)*sigPower/snr;

% Equalization 
eqSym = zeros(size(rxOFDM,1),numOFDMSymbols,numTx);
for i = 1:numOFDMSymbols
    eqSym(:,i,:) = ofdmEqualize(rxOFDM(:,i,:), ...
        squeeze(chanEst(:,i,:,:)),'Algorithm','ZF');
end

% Soft QAM demodulation
demodOut = qamdemod(eqSym, modOrder, 'UnitAveragePower', true, ...
	'OutputType', 'approxllr', 'NoiseVariance', nVar);
demodOut = demodOut(:);

% LDPC decoding
msgLen = cfgLDPCDec.NumInformationBits;
codewordLen = cfgLDPCDec.BlockLength;
numCodewordsPerFrame = length(demodOut)/codewordLen;
decBits = zeros(numCodewordsPerFrame*msgLen, 1);
for i = 1:numCodewordsPerFrame
    decBits((i-1)*msgLen+1:i*msgLen,:) = ...
        ldpcDecode(demodOut((i-1)*codewordLen+1:i*codewordLen,:),...
            cfgLDPCDec, 50, 'Termination', 'early');
end

end

% [EOF]