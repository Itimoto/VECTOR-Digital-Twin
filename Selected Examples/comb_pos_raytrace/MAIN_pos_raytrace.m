%% MIMO-OFDM Raytracing Comm + 802.11az Super-Resolution Time of Arrival Estimator
% Combination of MIMO-OFDM example: https://www.mathworks.com/help/comm/ug/indoor-mimo-ofdm-communication-link-using-ray-tracing.html
% and 802.11az Positioning example: https://www.mathworks.com/help/wlan/ug/802-11az-indoor-positioning-using-super-resolution-time-of-arrival-estimation.html
clear all;
close all;
%clc;

%% 3D Indoor Scenario
mapFileName = "conferenceroom.stl";

c = physconst('lightspeed');
fc = 5.8e9; % Carrier frequency
lambda = c/fc;

numSTAant = 1;  % # of Transmit Antennas (on User Terminal)
numAPant = 4;   % # of Receive Antennas (on Base Station)
txArray = arrayConfig("Size",[1 numSTAant],"ElementSpacing",2*lambda);
rxArray = arrayConfig("Size",[1 numAPant],"ElementSpacing",lambda);

helperViewArray(txArray); % To visualize the TX APA
helperViewArray(rxArray); % For RX URA APA

% Specify a transmitter site close to the upper corner of the room, which
% can be a WiFi accesss point
tx = txsite("cartesian", ...
    "Antenna",txArray, ...
    "AntennaPosition",[-1.46; -1.42; 2.1], ...
    'TransmitterFrequency',fc);
    %"AntennaPosition",[0; 25; 300],...% Empty space example

rx = rxsite("cartesian", ...
    "Antenna",rxArray, ...
    "AntennaPosition",[.3; .3; .85], ...
    "AntennaAngle",[0;90]);
    %"AntennaPosition",[0; 0; 0],...% Empty space example - set at origin

% Use the siteviewer function with the map file specified to view the scene
% in 3D. Use show function to visualize the transmitters and receivers
siteviewer("SceneModel",mapFileName);
%siteviewer("Terrain", 'none'); % To have 'nothing' / empty space
show(tx,"ShowAntennaHeight",false)
show(rx,"ShowAntennaHeight",false)

%% Ray Tracing Model
% Perform raytracing analysis between the transmitter and receiver sites
% and return the comm.Ray objects, using the Shooting and Bouncing rays
% (SBR) method. Specify the surface material of the scene as wood and
% search for rays with up to 2 reflections. SBR method supports up to 10
% orders of reflections
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",2, ...          % No reflections - LOS only
    "SurfaceMaterial","plasterboard");

rays = raytrace(tx,rx,pm);

% Extract computed rays from the cell array return
rays = rays{1,1};
% Examine the ray tracing results by looking at the number of reflections:
[rays.NumInteractions]
% ...propagation distance...
[rays.PropagationDistance]
% ...and path loss value of each ray
[rays.PathLoss]
% Use the plot function to plot the rays in the 3D scene in Siteviewer,
% with each ray colored based on path loss value
plot(rays,"Colormap",jet,"ColorLimits",[50, 95]);

%% Deterministic Channel Model from Ray Tracing
% Create a deterministic multipath channel model using the above ray
% tracing results. Specify the instananeous velocity of the receiver to
% reflect typical low mobility of a device in the indoor environment
rtChan = comm.RayTracingChannel(rays,tx,rx);
rtChan.SampleRate = 300e6;
rtChan.ReceiverVirtualVelocity = [0.1;0.1;0];

% Assign Eb/No value and derive SNR value from it for AWGN
bitsPerCarrier = 6; % Suppose we're using 64-QAM, which exists for 802.11ac & az
codeRate = 2/3;     % worst case 1/2, best case 5/6
EbNo = 30; % In dB
SNR = convertSNR(EbNo,"ebno", ...
  "BitsPerSymbol",bitsPerCarrier, ... % worst case 1, best case 10
  "CodingRate",codeRate);             % worst case 1/2, best case 5/6  
SNRLin = 10^(SNR/10);      % Linear

%% 802.11az Waveform Configuration
chanBW = "CBW20"; % 20MHz Channel
numSTS = min([numSTAant numAPant]);%numSTAant * numAPant; % # of Space-Time Streams
numLTFRepetitions = 8;  % # of HE-LTF repetitions

% Configure the HE ranging NDP parameters of the STA (User Terminal)
cfgSTABase = heRangingConfig;
cfgSTABase.ChannelBandwidth = chanBW;
cfgSTABase.NumTransmitAntennas = numSTAant;
cfgSTABase.SecureHELTF = true;
cfgSTABase.User{1}.NumSpaceTimeStreams = numSTS;
cfgSTABase.User{1}.NumHELTFRepetition = numLTFRepetitions;
cfgSTABase.GuardInterval = 1.6;

% Configure HE ranging NDP parameters of the AP (Base Station)
cfgAPBase = heRangingConfig;
cfgAPBase.ChannelBandwidth = chanBW;
cfgAPBase.NumTransmitAntennas = min([numSTAant numAPant]);
cfgAPBase.SecureHELTF = true;
cfgAPBase.User{1}.NumSpaceTimeStreams = numSTS;
cfgAPBase.User{1}.NumHELTFRepetition = numLTFRepetitions;
cfgAPBase.GuardInterval = 1.6;

ofdmInfo = wlanHEOFDMInfo('HE-LTF',chanBW,cfgSTABase.GuardInterval);
sampleRate = wlanSampleRate(chanBW);

% Miscellaneous channel config variables
chBaseInfo = info(rtChan);
chDelay = chBaseInfo.ChannelFilterDelay;
numPaths = size(rays, 2); % Number of paths simulated

%% Ranging Measurement
delayULDL = 16e-6; % Time Delay between UL NDP ToA and DL NDP ToD, in seconds

% Use a separate channel and waveform config object for each (potential)
% parfor stream, which would start here. -- Iterates over the SNR values
chan = rtChan; % Set it to the Raytracing Channel we generated earlier
cfgAP = cfgAPBase;
cfgSTA = cfgSTABase;

% Initialize ranging error and total failed packet count variables
rangingError = 0;
failedPackets = 0;

% Set random substream index per (potential) iteration to ensure each
% iteration uses a repeatable set of random numbers
stream = RandStream('combRecursive','Seed',654321);
stream.Substream = 1; % Would be equal to index in parfor loop
RandStream.setGlobalStream(stream);

% Define SNR per active subcarrier to account for noise energy in nulls
snrVal = SNR - 10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);

% 802.11az Pos Example starts Monte Carlo here
% Range-based delay
delay = distance(tx, rx)/c; % Divide out the speed of light
sampleDelay = delay*sampleRate;

% Loop over different APs starts here
linkType = ["Uplink","Downlink"];
% ToD of UL NDP (t1)
todUL = randsrc(1,1,0:1e-9:1e-6);
% Loop for both UL and DL transmission
numLinks = numel(linkType);
txTime = zeros(1,numLinks);

for l = 1:numLinks
    if linkType(l) == "Uplink" % STA to AP
        cfgSTA.UplinkIndication = 1; % For UL
        % Generate a random secure HE-LTF sequence for the exchange
        cfgSTA.User{1}.SecureHELTFSequence = dec2hex(randsrc(1,10,(0:15)))';
        cfg = cfgSTA;
    else                        % AP to STA
        % Generate a random secure HE-LTF sequence for the exchange
        cfgAP.User{1}.SecureHELTFSequence = dec2hex(randsrc(1,10,(0:15)))';
        cfg = cfgAP;    % For DL
    end

    % Set different channel for UL and DL, assuming that the channel is not
    % reciprocal
    reset(chan);

    % Generate HE Ranging NDP transmission
    txWaveform = heRangingWaveformGenerator(cfg);

    % Introduce time delay (fractional and integer) in the transmit waveform
    txDelay = heDelaySignal(txWaveform, sampleDelay);

    % Pad signal and pass through multipath channel (that's the raytraced channel!)
    [txMultipath, CIR] = chan(txDelay); % <-- This is the key change to the 802.11az example!

    % Pass waveform through AWGN channel
    rxWaveform = awgn(txMultipath,snrVal);

    % Perform synchronization and channel estimation
    [chanEstActiveSC,integerOffset] = heRangingSynchronize(rxWaveform,cfg); % << chanEstActiveSC looks like it could be CSI!

    % Estimate the transmission time between UL and DL
    if ~isempty(chanEstActiveSC) % If packet detection is successful
        % Estimate fractional delay with MUSIC super-resolution
        fracDelay = heRangingTOAEstimate(chanEstActiveSC,ofdmInfo.ActiveFFTIndices,...
                                        ofdmInfo.FFTLength, sampleRate, numPaths);
        integerOffset = integerOffset - chDelay;    % Account for channel filter delay
        intDelay = integerOffset/sampleRate;        % Estimate integer time delay
        txTime(l) = intDelay + fracDelay;           % Transmission time
    else % Packet detection failed
        txTime(l) = NaN;
    end
end

if ~any(isnan(txTime)) % If packet detection succeeds
    % TOA of UL waveform (t2)
    toaUL = todUL + txTime(1);

    % Time of departure of DL waveform (t3)
    todDL = toaUL + delayULDL;

    % TOA DL waveform (t4)
    toaDL = todDL + txTime(2);

    % Compute the RTT
    rtt = (toaDL-todUL) - (todDL-toaUL);

    % Estimate distance between STA and AP
    distEst = (rtt/2) * c;
    % Accumulate error
    rangingError = rangingError + abs(distance(tx, rx) - distEst);
else % If packet detection fails
    distEst = NaN;
    failedPackets = failedPackets + 1;
end

%% Display results
disp(['At SNR = ', num2str(SNR), ' dB, estimated distance ', ...
      num2str(distEst), 'm for true distance ', num2str(distance(tx, rx))]);

% Set up plots
figure; rows = 2; cols = 2; sgtitle("CSI Estimation from Ranging Process & Perfect Channel Estimator");
%% Plot CSI from generated by Ranging Process
ranging_chanEst = chanEstActiveSC; % Generated by heRangingSynchronize(...)
Hest = ranging_chanEst; Hest_flat = reshape(Hest, [size(Hest,1), size(Hest, 2)*size(Hest,3)]);
x = (-size(Hest_flat, 1)/2:1:(size(Hest_flat,1)/2 - 1))';
subplot(rows, cols, cols*0 + 1); plot(x, abs(Hest_flat)); title("Ranging: CSI Magnitude"); xlabel("Subcarrier Index?"); ylabel("Magnitude (dB?)"); grid on;
subplot(rows, cols, cols*1 + 1); plot(x, rad2deg(angle(Hest_flat))); title("Ranging: CSI Phase"); xlabel("Subcarrier Index?"); ylabel("Phase (deg)"); grid on; yticks(-180:30:180);

%% Plot CSI (?) via Perfect Channel Estimation
% From `helperIndoorRayTracingRxProcessing`, modified to fit the output
%  from ofdmInfo ~ wlanHEOFDMInfo (instead of direction OFDM mod/demodulator)
% Retrieve OFDM Parameters
fftLen = ofdmInfo.FFTLength;
cpLen = ofdmInfo.CPLength;
pilotCarrierIdx = ofdmInfo.PilotIndices;
%dataCarrierIdx = setdiff( ... %<<<--- Original
%    numGuardBandCarriers(1)+1:fftLen-numGuardBandCarriers(2), ...
%    [pilotCarrierIdx; fftLen/2+1]);
dataCarrierIdx = ofdmInfo.DataIndices; % Use it directly for data subcarriers
% Data carriers excluding pilots and guard bands, so no need for setdiff
% However, if DC carrier needs to be excluded, keep the following line:
dataCarrierIdx = setdiff(dataCarrierIdx, fftLen/2 + 1);

% Perfect Channel Estimation
chanDelay = channelDelay(CIR, chBaseInfo.ChannelFilterCoefficients);
chanEst = ofdmChannelResponse(CIR, chBaseInfo.ChannelFilterCoefficients, ...
                                fftLen, cpLen, dataCarrierIdx, chanDelay);
% Plot.
% Case when exposing chanEst, which has a bunch of symbols. We pick the first symbol, then graph the CSI for it
dims = size(chanEst(:,1,:,:)); Hest = reshape(chanEst(:,1,:,:), [dims(1),dims(3:end)]);
Hest_flat = reshape(Hest, [size(Hest,1), size(Hest, 2)*size(Hest,3)]);
x = (-size(Hest_flat, 1)/2:1:(size(Hest_flat,1)/2 - 1))';
subplot(rows, cols, cols*0 + 2); plot(x, abs(Hest_flat)); title("Perfect: CSI Magnitude"); xlabel("Subcarrier Index?"); ylabel("Magnitude (dB?)"); grid on;
subplot(rows, cols, cols*1 + 2); plot(x, rad2deg(angle(Hest_flat))); title("Perfect: CSI Phase"); xlabel("Subcarrier Index?"); ylabel("Phase (deg)"); grid on; yticks(-180:30:180);