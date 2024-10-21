function heRangingPlotDelayProfile
%heRangingPlotDelayProfile Plots the estimated multipath delay profile
%
%   heRangingPlotDelayProfile simulates a single 802.11az packet
%   transmission through a TGax multipath channel and estimates the delay
%   profile and the time of arrival (ToA) via MUSIC super-resolution.

%   Copyright 2020-2021 The MathWorks, Inc.

% Configure simulation parameters
snr = 35; % SNR (dB)
distance = 10; % Transmit - receive distance in meters

s = rng(); % Store current random state
rng(99,'combRecursive'); % Set repeatable random number stream

% Configure the HE ranging NDP parameters of the STA
chanBW = 'CBW40'; % Channel bandwidth
numTx = 1; % Number of transmit antennas
numRx = 1; % Number of receive antennas

cfgSTA = heRangingConfig;
cfgSTA.ChannelBandwidth = chanBW;
cfgSTA.NumTransmitAntennas = numTx;
cfgSTA.SecureHELTF = true;
cfgSTA.User{1}.NumSpaceTimeStreams = 1;
cfgSTA.User{1}.NumHELTFRepetition = 3;
cfgSTA.GuardInterval = 1.6;
cfgSTA.UplinkIndication = 1; % For UL
cfgSTA.User{1}.SecureHELTFSequence = dec2hex(randsrc(1,10,(0:15)))';

ofdmInfo = wlanHEOFDMInfo('HE-LTF',chanBW,cfgSTA.GuardInterval);
sampleRate = wlanSampleRate(chanBW);
speedOfLight = physconst('lightspeed'); % m/s

% Configure multipath channel
chan = wlanTGaxChannel;
chan.DelayProfile = 'Model-B';
chan.NumTransmitAntennas = numTx;
chan.NumReceiveAntennas = numRx;
chan.SampleRate = sampleRate;
chan.CarrierFrequency = 5e9;
chan.ChannelBandwidth = chanBW;
chan.PathGainsOutputPort = true;
chan.NormalizeChannelOutputs = false;

chBaseInfo = info(chan);
chDelay = chBaseInfo.ChannelFilterDelay;
numPaths = size(chBaseInfo.PathDelays,2);

% Generate 11az waveform
tx = heRangingWaveformGenerator(cfgSTA);

% Introduce time delay (fractional and integer) in the transmit waveform
delay = distance/speedOfLight;
sampleDelay = sampleRate*delay;
txDelay = heDelaySignal(tx,sampleDelay);

% Pad signal and pass through multipath channel
[txMultipath,pathGains] = chan([txDelay;zeros(50,cfgSTA.NumTransmitAntennas)]);

% Pass the waveform through AWGN channel
% Account for noise energy in nulls so the SNR is defined per
% active subcarrier
rx = awgn(txMultipath,snr - 10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones));

% Perform synchronization and channel estimation
[chanEstActiveSC,integerOffset] = heRangingSynchronize(rx,cfgSTA);

% Estimate fractional delay w. MUSIC super-resolution
[fracDelay,delProfTimeEst,delProfPowEst] = heRangingTOAEstimate(chanEstActiveSC, ...
                                                                  ofdmInfo.ActiveFFTIndices, ...
                                                                  ofdmInfo.FFTLength,...
                                                                  sampleRate,...
                                                                  numPaths);

% Account for channel filter delay
integerOffset = integerOffset-chDelay;

% Plot the estimated delay profile
truePathDelays = delay+chBaseInfo.PathDelays;
truePathPowers = abs(pathGains(1,:,1,1)).^2;
estPathDelays = delProfTimeEst+integerOffset/sampleRate;
estToa = fracDelay+integerOffset/sampleRate;

figure;
stem(1e9*truePathDelays,truePathPowers);
K = truePathPowers(1)/delProfPowEst(estPathDelays==estToa);
hold on,plot(1e9*estPathDelays,K*delProfPowEst);
hold on,plot(1e9*estToa,K*delProfPowEst(estPathDelays==estToa),'r*');
xlim([0 1.1e9*truePathDelays(end)])
legend('True delay profile','MUSIC estimate','Estimated ToA')
xlabel('Time (ns)'),ylabel('Path power')
title('Delay profile (40 MHz, Model-B, 35 dB SNR, SISO)')

% Restore random state
rng(s);

end
