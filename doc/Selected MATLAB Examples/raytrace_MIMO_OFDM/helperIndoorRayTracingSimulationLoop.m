function ber = helperIndoorRayTracingSimulationLoop(cfgLDPCEnc, cfgLDPCDec, ofdmMod, ofdmDemod, ...
    rtChan, errRate, modOrder, numCodewordsPerFrame, EbNo)
%HELPERINDOORRAYTRACINGSIMULATIONLOOP Calculate BER for a given set of Eb/No values

%   Copyright 2020-2021 The MathWorks, Inc. 

% Simulation setup
maxNumFrms = 50;
maxNumErrs = 300;

% Retrieve parameters
codeRate = cfgLDPCEnc.CodeRate;
bitsPerSymbol = log2(modOrder);
SNR = convertSNR(EbNo,"ebno", ...
  "BitsPerSymbol",bitsPerSymbol, ...
  "CodingRate",codeRate);
SNRLin = 10.^(SNR/10);      % Linear
sigPower = ofdmMod.NumTransmitAntennas/ofdmDemod.NumReceiveAntennas;
rtChanInfo = info(rtChan);

% Set RNG for repeatability
rng(100); 

% Initialize BER
ber = zeros(3,length(EbNo));

% Simulation loop
for i = 1:length(EbNo)
    reset(errRate);
    reset(rtChan);
    frmIdx = 0; 
    while ((ber(2,i) < maxNumErrs) && (frmIdx < maxNumFrms))
        % Tx processing
        [txWave, srcBits] = helperIndoorRayTracingWaveformGen( ...
            numCodewordsPerFrame, cfgLDPCEnc, modOrder, ofdmMod);
        % Channel
        chanIn = [txWave; ....
            zeros(ofdmMod.FFTLength + ofdmMod.CyclicPrefixLength, ...
            ofdmMod.NumTransmitAntennas)];
        [chanOut, CIR] = rtChan(chanIn);
        rxWave = awgn(chanOut, SNRLin(i), sigPower, 'linear');
        % Rx processing
        decBits = helperIndoorRayTracingRxProcessing( ...
            rxWave, CIR, rtChanInfo, cfgLDPCDec, modOrder, ofdmDemod, SNRLin(i));
        % Update BER
        ber(:,i) = errRate(srcBits, double(decBits));
        
        frmIdx = frmIdx + 1;
    end
end

% Plot BER vs. EbNo figure
figure;
semilogy(EbNo, ber(1,:), 'o-');
grid on;
xlabel('Eb/No (dB)'); ylabel('BER');

end

% [EOF]