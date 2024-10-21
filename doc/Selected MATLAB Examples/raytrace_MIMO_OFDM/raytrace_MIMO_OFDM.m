%% Indoor MIMO-OFDM Communication Link using Ray Tracing
% https://www.mathworks.com/help/comm/ug/indoor-mimo-ofdm-communication-link-using-ray-tracing.html
% Unlike stocahstic models, the ray-tracing method is 3-D environment and
% transceiver sites specific, and can have high sensitivity in the
% surrounding environment.

%% 3D Indoor Scenario
mapFileName = "conferenceroom.stl";

fc = 5.8e9; % Carrier frequency
lambda = physconst("lightspeed")/fc;

% The TX antenna is a 4-element ULA with twice the wavelength between
% elements. RX is a 4-element URA w/ 1 wavelength between elements
txArray = arrayConfig("Size",[4,1],"ElementSpacing",2*lambda);
rxArray = arrayConfig("Size",[4 4],"ElementSpacing",lambda);

helperViewArray(txArray); % To visualize the TX APA
helperViewArray(rxArray); % For RX URA APA

% Specify a transmitter site close to the upper corner of the room, which
% can be a WiFi accesss point
tx = txsite("cartesian", ...
    "Antenna",txArray, ...
    "AntennaPosition",[-1.46; -1.42; 2.1], ...
    'TransmitterFrequency',5.8e9);

rx = rxsite("cartesian", ...
    "Antenna",rxArray, ...
    "AntennaPosition",[.3; .3; .85], ...
    "AntennaAngle",[0;90]);

% Use the siteviewer function with the map file specified to view the scene
% in 3D. Use show function to visualize the transmitters and receivers
siteviewer("SceneModel",mapFileName);
show(tx,"ShowAntennaHeight",false)
show(rx,"ShowAntennaHeight",false)

%% Ray Tracing
% Perform raytracing analysis between the transmitter and receiver sites
% and return the comm.Ray objects, using the Shooting and Bouncing rays
% (SBR) method. Specify the surface material of the scene as wood and
% search for rays with up to 2 reflections. SBR method supports up to 10
% orders of reflections
pm = propagationModel("raytracing", ...
    "CoordinateSystem","cartesian", ...
    "Method","sbr", ...
    "AngularSeparation","low", ...
    "MaxNumReflections",2, ...
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
rtChan.ReceiverVirtualVelocity = [0.1;0.1;0]

% Use the `showProfile` object function to visualize the:
% Power Delay Profile (PDP), Angle of Departure (AoD), and 
%   Angle of Arrival (AoA) of the rays in the channel
% (the visualization has PDP taking into account the transmit/receive array
% pattern gains in addition to the path loss for each ray)
showProfile(rtChan);
% Use the `info` object function to obtain the number of transmit and
% receive elements
rtChanInfo = info(rtChan)

numTx = rtChanInfo.NumTransmitElements;
numRx = rtChanInfo.NumReceiveElements;

%% System Parameters
% Configure a comms link using LDPC coding, 64 QAM and OFDM with 256
% subcarriers
% Specify 4 LDPC codewords per frame, which results in 50 OFDM symbols per frame
% Create LDPC encoder and decoder configuration objects
cfgLDPCEnc = ldpcEncoderConfig(dvbs2ldpc(1/2));
cfgLDPCDec = ldpcDecoderConfig(cfgLDPCEnc);
numCodewordsPerFrame = 4;
codewordLen = cfgLDPCEnc.BlockLength;

% Parameters for QAM modulation per subcarrier
bitsPerCarrier = 6;
modOrder = 2^bitsPerCarrier;
codeRate = cfgLDPCEnc.CodeRate;

% Create OFDM modulator and demodulator objects 
fftLen = 256; 
cpLen = fftLen/4; 
numGuardBandCarriers = [9; 8];
pilotCarrierIdx = [19:10:119, 139:10:239]';
numDataCarriers = ...
    fftLen - sum(numGuardBandCarriers) - length(pilotCarrierIdx) - 1;
numOFDMSymbols = ...
    numCodewordsPerFrame * codewordLen / ...
    bitsPerCarrier / numDataCarriers / numTx;
ofdmMod = comm.OFDMModulator( ...
    "FFTLength",fftLen, ....
    "NumGuardBandCarriers",numGuardBandCarriers, ...
    "InsertDCNull",true, ...
    "PilotInputPort",true, ...
    "PilotCarrierIndices",pilotCarrierIdx, ...
    "CyclicPrefixLength",cpLen, ...
    "NumSymbols",numOFDMSymbols, ...
    "NumTransmitAntennas",numTx);
ofdmDemod = comm.OFDMDemodulator(ofdmMod);
ofdmDemod.NumReceiveAntennas = numRx;
cd = comm.ConstellationDiagram( ...    
    "ReferenceConstellation", qammod(0:modOrder-1, modOrder, 'UnitAveragePower', true), ...
    "XLimits", [-2 2], ...
    "YLimits", [-2 2]);

% Error rate calculation object to compute bit error rate (BER)
errRate = comm.ErrorRate;
% Assign Eb/No value and derive SNR value from it for AWGN
EbNo = 30; % In dB
SNR = convertSNR(EbNo,"ebno", ...
  "BitsPerSymbol",bitsPerCarrier, ...
  "CodingRate",codeRate);
SNRLin = 10^(SNR/10);      % Linear

%% Link Simulation
% `helperIndoorRayTracingWaveformGen` generates a waveform consisting of
% one frame at the transmitter site by performing the following:
%   1. Encode randomly generated bits by LDPC
%   2. Modulate encoded bits by 64-QAM
%   3. Apply OFDM modulation to convert signals from frequency domain to
%       time domain
rng(100); % Set RNG for repeatability
[txWave,srcBits] = ...
    helperIndoorRayTracingWaveformGen( ...
    numCodewordsPerFrame,cfgLDPCEnc,modOrder,ofdmMod);
% Pass waveform through the Ray Tracing Channel Model and add white noise
% To account for channel filtering delay, append an additional null OFDM
% symbol to the end of the waveform
chanIn = [txWave; zeros(fftLen + cpLen, numTx)];
[chanOut,CIR] = rtChan(chanIn);
rxWave = awgn(chanOut,SNRLin,numTx/numRx,'linear');

% `helperIndoorRayTracingRxProcessing` function decodes the
% channel-impaired waveform at the receiver by performing:
%`  1. Perfect channel estimation using the Channel Impulse Response (CIR)
%       output and the channel filter coefficients from the channel object's info
%       method
%   2. OFDM demodulation to bring the signals back to the frequency domain
%   3. Symbol equalization on each subcarrier
%   4. Soft 64-QAM demodulation to get LLR
%   5. LDPC Decoding
[decBits, eqSym] = ...
    helperIndoorRayTracingRxProcessing(rxWave,CIR, ...
    rtChanInfo,cfgLDPCDec,modOrder,ofdmDemod,SNRLin);
cd(eqSym(:));

% Find BER
ber = errRate(srcBits,double(decBits));
disp(ber(1));

% To plot BER curve against of range of EbNo values...:
%EbNoRange = 27:36;
%helperIndoorRayTracingSimulationLoop( ...
%    cfgLDPCEnc,cfgLDPCDec,ofdmMod,ofdmDemod,rtChan,errRate, ...
%    modOrder,numCodewordsPerFrame,EbNoRange);