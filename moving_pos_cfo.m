clear all;
%% DOA ESTIMATOR - MOVING POSITION - CARRIER FREQUENCY OFFSET APPLIED
CFO = 1000; % (Hz) - TX to RX
PERM_ARRAY = [[2 3];[1 4]];
INTER_CFO = [-1000, 2500];
%INTER_CFO = [0, 0];
MIN_SWITCH_TIME = inf; % Flip a coin every this many seconds. `inf` to turn off.

EBNO_IN = 10;
NUM_REFLECTIONS = 10; % 0 for LoS

%% MIMO-OFDM Raytracing Comm CSI Generator
% Notes - Specify 
% Combination of MIMO-OFDM example: https://www.mathworks.com/help/comm/ug/indoor-mimo-ofdm-communication-link-using-ray-tracing.html
% and 802.11az Positioning example: https://www.mathworks.com/help/wlan/ug/802-11az-indoor-positioning-using-super-resolution-time-of-arrival-estimation.html
%close all;
%clc;

%% Array Configuration
c = physconst('lightspeed');
fc = 2.447e9; % Carrier Frequency
lambda = c/fc;
spacing = 0.079961058/lambda;

numSTAant = 1; % # of Transmit Antennas (on User Terminal / STA)
numAPant  = 4; % # of Receive Antennas (on Base Station / AP)
txArray = arrayConfig("Size",[1 numSTAant],"ElementSpacing",spacing*lambda);
rxArray = arrayConfig("Size",[1 numAPant],"ElementSpacing",spacing*lambda);

% Use the siteviewer function with the map file specified to view the scene
% in 3D. Use show function to visualize the transmitters and receivers
%filename = "env-models/Copernicuscrater3Xv.stl"; % (pulled from https://nasa3d.arc.nasa.gov/detail/copernicus-crater)
%siteviewer(SceneModel=filename, ShowEdges=false); COORDSYS="cartesian"; % Uncomment for custom STL
%siteviewer("Terrain", 'none', "Hidden",true); COORDSYS="cartesian";% To have 'nothing' / empty space
siteviewer("SceneModel","conferenceroom.stl"); COORDSYS="cartesian";% Uncomment to simulate multipath
%site = siteviewer("Hidden",true, Basemap="openstreetmap", Buildings="manhattan.osm"); COORDSYS="geographic";

% Specify a transmitter site close to the upper corner of the room, which
% can be a WiFi accesss point
tx = txsite(...
    COORDSYS, ...%"geographic" or "cartesian"
    "Antenna",txArray, ...
    'TransmitterFrequency',fc,...
    'AntennaHeight',1);

rx = rxsite( ...
    COORDSYS,...
    "Antenna",rxArray, ...
    "AntennaAngle",[0; 0],...
    "AntennaHeight",1);

%% Physical Locations
% R_a    ~ Linear Distance from TX to RX (m)
% beta_a ~ Direction of TX, with respect to RX (in degrees)
%           0 degrees is parallel to the array, on the right-hand side
%               relative to the X-Z plane visualized by `helperViewArray`
%           10 degrees moves clockwise about the Z axis
% dt     ~ Time step in between each snapshot (s)
ft2m = 0.3048; % Feet to Meters
startPt = 5; endPt = 3;
R_a = [repmat([startPt*ft2m], 1, 50)];
b_a = repmat([72], 1, length(R_a)); % Just have it run at boresight for now
%b_a = [90 90 90 91 92 93 94 95 96 97 98 99 100 100 100 100 100 100 100];
%R_a = repmat([1], 1, length(b_a));
dt  = 0.23;
timestamps = (1:length(R_a))*dt; % For tracking/inducing Carrier Frequency Offset

%% Additional Channel Parameters (e.g. SNR)
% Assign Eb/No value and derive SNR value from it for AWGN
bitsPerCarrier = 6; % Suppose we're using 64-QAM, which exists for 802.11ac & az
codeRate = 2/3;     % worst case 1/2, best case 5/6
EbNo = EBNO_IN; % In dB
SNR = convertSNR(EbNo,"ebno", ...
  "BitsPerSymbol",bitsPerCarrier, ... % worst case 1, best case 10
  "CodingRate",codeRate);             % worst case 1/2, best case 5/6  
SNRLin = 10^(SNR/10);      % Linear

%% 802.11az Waveform Configuration
bw = 20;
chanBW = "CBW"+num2str(bw);
numSTS = min([numSTAant numAPant]); % # of Space-Time Streams
numLTFRepetitions = 8; % # of HE-LTF repetitions

% Configure the HE ranging NDP parameters of the STA (User Terminal)
cfgSTABase = ranging.heRangingConfig;
cfgSTABase.ChannelBandwidth = chanBW;
cfgSTABase.NumTransmitAntennas = numSTAant;
cfgSTABase.SecureHELTF = true;
cfgSTABase.User{1}.NumSpaceTimeStreams = numSTS;
cfgSTABase.User{1}.NumHELTFRepetition = numLTFRepetitions;
cfgSTABase.GuardInterval = 0.8;

% Configure HE ranging NDP parameters of the AP (Base Station)
cfgAPBase = ranging.heRangingConfig;
cfgAPBase.ChannelBandwidth = chanBW;
cfgAPBase.NumTransmitAntennas = min([numSTAant numAPant]);
cfgAPBase.SecureHELTF = true;
cfgAPBase.User{1}.NumSpaceTimeStreams = numSTS;
cfgAPBase.User{1}.NumHELTFRepetition = numLTFRepetitions;
cfgAPBase.GuardInterval = 0.8;

ofdmInfo = wlanHEOFDMInfo('HE-LTF',chanBW,cfgSTABase.GuardInterval);
sampleRate = wlanSampleRate(chanBW);
subcarrierSpacing = sampleRate/ofdmInfo.FFTLength;

%% Ranging Setup
delayULDL = 16e-6; % Time Delay between UL NDP ToA and DL NDP ToD, in seconds

cfgAP = cfgAPBase;
cfgSTA = cfgSTABase;

% Initialize ranging error and total failed packet count variables
rangingError = 0;
failedPackets = 0;

% Set random substream index per (potential) iteration to ensure each
% iteration sues a repeatable set of random numbers
stream = RandStream('combRecursive','Seed',654321);
stream.Substream = 1; % Equal to index in parfor loop
RandStream.setGlobalStream(stream);

% Define SNR per active subcarrier to account for noise energy in nulls
snrVal = SNR - 10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);

linkType = ["Uplink","Downlink"];
% ToD of UL NDP (t1)
todUL = randsrc(1,1,0:1e-9:1e-6);
% Loop for both UL and DL transmission
numLinks = numel(linkType);
txTime = zeros(1,numLinks);

numSnapshots = length(R_a);
Hest_Multi = zeros(numSTAant, numAPant, ofdmInfo.NumTones, numSnapshots); % [AT AR S K]
subcFreq_Multi = zeros(ofdmInfo.NumTones, numSnapshots); % [S K] (In case it changes?)
disp("Starting simulation...")
for k = 1:numSnapshots
    disp("Simulating Snapshot "+string(k)+" of "+string(numSnapshots));
    %% Place TX/RX Down:
    if COORDSYS == "cartesian"
        rx.AntennaPosition = [0; 0; 0.3]; % Place RX Antenna at the origin   
        %rx.AntennaPosition = [287.5; 0.2; 5.305]; % Place RX Antenna at the origin of the lunar STL   
        tx.AntennaPosition = [R_a(k)*cos(deg2rad(b_a(k) - 90)); R_a(k)*sin(deg2rad(b_a(k) - 90)); -(0.2-0.3)];%-0.005];% Place TX Antenna some distance away
        tx.AntennaPosition = tx.AntennaPosition + rx.AntennaPosition; % Place TX relative to rx
    elseif COORDSYS == "geographic"
        rx.Latitude = 40.7089; rx.Longitude = -74.009; rx.AntennaAngle = [0 0];
        [tx.Latitude, tx.Longitude] = helper.offsetLatLon(rx.Latitude, rx.Longitude, R_a(k), 180-b_a(k)); % Assuming RX is pointed northward.
    end

    % Range-based delay
    delay = distance(tx, rx)/c; % Divide out the speed of light
    sampleDelay = delay*sampleRate;

    %% Set up Raytracing Model (Shooting & Bouncing Rays / SBR)
    pm = propagationModel("raytracing", ...
    "CoordinateSystem",COORDSYS, ...
    "Method","sbr", ...
    "AngularSeparation","low", ...
    "SurfaceMaterial","concrete");
    pm.MaxNumReflections = NUM_REFLECTIONS; % If 0: No reflections, LOS only
    
    rays = raytrace(tx,rx,pm);
    
    % Extract computed rays from the cell array return
    rays = rays{1,1};

    %% Display for User:
    %{
    show(tx, "ShowAntennaHeight",false)
    show(rx, "ShowAntennaHeight",false)
    %pattern(rx, fc); 
    plot(rays,"Colormap",jet,"ColorLimits",[50, 95]);
    %input("Good?");
    %}

    %% Deterministic Channel Model
    rtChan = comm.RayTracingChannel(rays,tx,rx);
    rtChan.SampleRate = 300e6;
    
    % Determine Receiver Virtual Velocity (for Doppler effects)
    if k == 1 % Start at rest
        rtChan.ReceiverVirtualVelocity = [0;0;0];
    else      % Use the movement from the last point to the current point
        ds = [R_a(k) - R_a(k-1)             ;   b_a(k) - b_a(k-1)             ; 0];
        ds = [ds(1)*cos(deg2rad(ds(2) - 90));   ds(1)*sin(deg2rad(ds(2) - 90)); 0]; % Convert to Cartesian
        rtChan.ReceiverVirtualVelocity = ds/dt; % change in position over time
    end

    % Miscellaneous channel config variables
    chBaseInfo = info(rtChan);
    chDelay = chBaseInfo.ChannelFilterDelay;
    numPaths = size(rays, 2); % Number of paths simulated
    chan = rtChan;

    %% Simulate the Link:
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
        txWaveform = ranging.heRangingWaveformGenerator(cfg);
    
        % Introduce time delay (fractional and integer) in the transmit waveform
        txDelay = ranging.heDelaySignal(txWaveform, sampleDelay);
    
        % Introduce Center Frequency Offset in the transmit waveform
        txCFO = nonideal.freqOffsetSignalWithTime(txDelay, CFO, sampleRate, (k-1)*dt);

        % Pad signal and pass through multipath channel (that's the raytraced channel!)
        [txMultipath, CIR] = chan(txCFO); % <-- This is the key change to the 802.11az example!

        % Introduce Inter-NIC CFO between the outer two traces (connected to the NIC)
        txInterCFO = nonideal.applyInterCFO(txMultipath, PERM_ARRAY, INTER_CFO, sampleRate, (k-1)*dt);

        % Swap 'Antennas' Internally
        txSwapped = nonideal.swapTracesRandom(txInterCFO, PERM_ARRAY, MIN_SWITCH_TIME, sampleRate, (k-1)*dt);

        % Pass waveform through AWGN channel
        rxWaveform = awgn(txSwapped,snrVal);
    
        % Perform synchronization and channel estimation
        [chanEstActiveSC,integerOffset] = ranging.heRangingSynchronize(rxWaveform,cfg); % << chanEstActiveSC looks like it could be CSI!

        % Output Per-Snapshot Channel State Information [AT AR S K]
        if (size(chanEstActiveSC, 1) > 0) % Sometimes, we can't solve for it. This frame is dropped.
            Hest_Multi(:, :, :, k) = permute(chanEstActiveSC, [2 3 1]); % [AT AR S K]
            subcFreq_Multi(:, k) = fc + ofdmInfo.ActiveFrequencyIndices*subcarrierSpacing;
        end
                
        % Estimate the transmission time between UL and DL
        if ~isempty(chanEstActiveSC) % If packet detection is successful
            % Estimate fractional delay with MUSIC super-resolution
            fracDelay = ranging.heRangingTOAEstimate(chanEstActiveSC,ofdmInfo.ActiveFFTIndices,...
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
end

%% Determine DOA via MUSIC Method
disp("Starting DOA Estimation...");
Hest = Hest_Multi;
Hest(Hest(1, 1, 1, :) == 0) = []; % Remove dropped frames (if any were found)
subcFreq = subcFreq_Multi(:, 1); % Use the first, assume no change occurs.
rxElemPos = helper.getElementPosition(rxArray).'; % Extract Element Positions
musicAvgWindow = 2;
thetaRange = [65 115];

[doa_music, doaMat_music, doaMUSIC] = doa_lib.naive_music(Hest, subcFreq, rxElemPos, musicAvgWindow, thetaRange);
% Plot DOA from MUSIC:
plotting.plotDoA_overSnapshots(doaMUSIC, "DOA Estimate via MUSIC", [65 115], musicAvgWindow);
plotting.plotDoA_overSubcarriers(doaMUSIC, "DOA Estimate via MUSIC", [65 115]);

% Set up plots
%% Plot CSI from generated by Ranging Process
plotting.plotCSI(Hest_Multi, "CSI Estimation from Ranging Process", subcFreq_Multi(:, 1));

%% Plot Rays in Siteviewer:
%siteviewer(Basemap="openstreetmap", Buildings="manhattan.osm");
show(tx, "ShowAntennaHeight",false)
show(rx, "ShowAntennaHeight",false)
pattern(rx, fc); 
plot(rays,"Colormap",jet,"ColorLimits",[50, 95]);

% Save CSI to simulate with VECTOR System
if upper(input("Save CSI? (y for yes): ", "s")) == 'Y'
    centerFreq = fc;                                    % Legacy
    chanBW = bw*1e6;                                    % Legacy
    elemPos = helper.getElementPosition(rxArray).';     % Element Position for DOA
    outputMatrix = Hest_Multi;                          % [AT AR S K]
    subcFreq = subcFreq_Multi(:, 1);                    % For DOA and Ranging
    timestamps = timestamps;                            % In seconds

    filename = input("Please enter new filename: ", 's');
    save(string(filename)+'.mat', 'centerFreq', 'chanBW', 'elemPos', 'outputMatrix', 'subcFreq', 'timestamps');
end

%{
diffs = zeros(numSnapshots-1, 1);
for k = 2:numSnapshots
    diffs(k) = angle(Hest_Multi(1, 1, 1, k)) - angle(Hest_Multi(1, 1, 1, k-1));
end
figure; plot(unwrap(diffs));
%}