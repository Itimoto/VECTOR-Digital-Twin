%% 802.11az Positioning Using Super-Resolution Time of Arrival Estimation
% https://www.mathworks.com/help/wlan/ug/802-11az-indoor-positioning-using-super-resolution-time-of-arrival-estimation.html
% Shows how to estimate the position of a station (STA) in a multipath
%   environment using a Time-of-Arrival (ToA) based positioning algorithm
%   defined in IEEE 802.11az WiFi standard
% Example estimates the ToA by using a multiple signal classification
%   (MUSIC) super-resolution approach, then estimates the two-dimensional
%   position of a STA by using trilateration
% Sim evaluates & compares the performance of the positioning algorithm at
%   multiple signal-to-noise ratio (SNR) points

% 802.11az supports two high-efficiency (HE) ranging Physical Layer (PHY)
%   Protocol Data Unit (PPDU) formats (analogous to the 802.11ax versions)
% - HE ranging null data packet (NDP)
% - HE trigger-based (TB) ranging NDP
% To parameterize & generate HE ranging NDPs, see the 
%   802.11az Waveform Generation Example

%{ 
    Modeling the measurement exchange btwn STA & APs:
    1. Generate a ranging NDP
    2. Delay the NDP according to a randomly generate distance between the
        STA & AP, adding fractional and integer sample delay
    3. Pass the waveform through an indoor TGax channel. 
        The example models different channel realizations for different
        packets
    4. Add Additive White Gaussian Noise (AWGN) to the received waveform
        The example uses the same SNR value for all links btwn STA & APs
    5. Perform synchronization and frequency correction on the received
        waveform
    6. Demodulate the HE-LTF (Long Training Field)
    7. Estimate the Channel Frequency Response from the HE-LTF
    8. Estimate the distance by using the MUSIC super-resolution algorithm
    9. Combine distance estimates from other STA-AP pairs and trilaterate
        pos'n of the STA
%}

%% Simulation Parameters
numIterations = 50;
snrRange = 15:10:35; % SNR points in dB
numAPs = 3; % Number of APs

%% 802.11az Waveform Configuration
% Waveform Generators for each AP & STA
chanBW = "CBW20";
numTx = 1; % # of Transmit Antennas
numRx = 4; % Number of Receive ANtennas
numSTS = 1; % Number of Space-Time Streams
numLTFRepetitions = 8; % Number of HE-LTF repetitions

% Configure the HE ranging NDP parameters of the STA
cfgSTABase = heRangingConfig;
cfgSTABase.ChannelBandwidth = chanBW;
cfgSTABase.NumTransmitAntennas = numTx;
cfgSTABase.SecureHELTF = true;
cfgSTABase.User{1}.NumSpaceTimeStreams = numSTS;
cfgSTABase.User{1}.NumHELTFRepetition = numLTFRepetitions;
cfgSTABase.GuardInterval = 1.6;

% Configure the HE ranging NDP parameters of the APs
cfgAPBase = cell(1,numAPs);
for iAP = 1:numAPs
    cfgAPBase{iAP} = heRangingConfig;
    cfgAPBase{iAP}.ChannelBandwidth = chanBW;
    cfgAPBase{iAP}.NumTransmitAntennas = numTx;
    cfgAPBase{iAP}.SecureHELTF = true;
    cfgAPBase{iAP}.User{1}.NumSpaceTimeStreams = numSTS;
    cfgAPBase{iAP}.User{1}.NumHELTFRepetition = numLTFRepetitions;
    cfgAPBase{iAP}.GuardInterval = 1.6;
end

ofdmInfo = wlanHEOFDMInfo('HE-LTF',chanBW,cfgSTABase.GuardInterval);
sampleRate = wlanSampleRate(chanBW);

%% Channel Configuration
delayProfile = "Model-B";

carrierFrequency = 5e9; % fc in Hz
speedOfLight = physconst('lightspeed');

chanBase = wlanTGaxChannel;
chanBase.DelayProfile = delayProfile;
chanBase.NumTransmitAntennas = numTx;
chanBase.NumReceiveAntennas = numRx;
chanBase.SampleRate = sampleRate;
chanBase.CarrierFrequency = carrierFrequency;
chanBase.ChannelBandwidth = chanBW;
chanBase.PathGainsOutputPort = true;
chanBase.NormalizeChannelOutputs = false;
% Get channel filter delay & number of paths
chBaseInfo = info(chanBase);
chDelay = chBaseInfo.ChannelFilterDelay;
numPaths = size(chBaseInfo.PathDelays, 2);

%% Ranging Measurement
delayULDL = 16e-6; % Time Delay between UL NDP ToA and DL NDP ToD, in seconds

numSNR = numel(snrRange);
distEst = zeros(numAPs, numIterations, numSNR); % Estimated Distance
distance = zeros(numAPs, numIterations, numSNR); % True Distance
positionSTA = zeros(2, numIterations, numSNR); % Two-dimensional position of the STA
positionAP = zeros(2, numAPs, numIterations, numSNR); % Two-dimensional positions of the APs
per = zeros(numSNR, 1); % Packet Error Rate (PER)

%parfor isnr = 1:numSNR % Use 'parfor' to speed up the simulation
for isnr = 1:numSNR
    
    % Use a separate channel and waveform configuration object for each parfor stream
    chan = chanBase;
    cfgAP = cfgAPBase;
    cfgSTA = cfgSTABase;
    
    % Initialize ranging error and total failed packet count variables
    rangingError = 0;
    failedPackets = 0;
    
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',654321);
    stream.Substream = isnr;
    RandStream.setGlobalStream(stream);

    % Define the SNR per active subcarrier to account for noise energy in nulls
    snrVal = snrRange(isnr) - 10*log10(ofdmInfo.FFTLength/ofdmInfo.NumTones);
    
    for iter = 1:numIterations
        
        % Position everyone:
        % Put STA away from AP
        posSTA = [0; 10]; % Positioning STA away from origin
        positionSTA(:,iter,isnr) = posSTA; 

        % Set up AP Positions
        arrayAperture = numAPs*6*speedOfLight/carrierFrequency; % Rel to Lambda size
        posAP = [linspace((1e-3 * rand())-arrayAperture/2, (1e-3 * rand())+arrayAperture/2, numAPs); 
                                     [(1e-3 * rand())+0*ones(1,3)]]; % X & Y positions
                                    % Y is at 0 for each (ULA)
                                    % X is spread evenly across array aperture
        positionAP(:,:,iter,isnr) = posAP;

        % Distance from each AP to the STA
        distanceAllAPs = sqrt((posAP(1,:) - posSTA(1)).^2 + (posAP(2,:) - posSTA(2)).^2);
        distance(:,iter,isnr) = distanceAllAPs;

        % OG Position Generation:
        %[positionSTA(:,iter,isnr),positionAP(:,:,iter,isnr),distanceAllAPs] = heGeneratePositions(numAPs);
        %distance(:,iter,isnr) = distanceAllAPs;
        
        % Range-based delay
        delay = distance(:,iter,isnr)/speedOfLight;
        sampleDelay = delay*sampleRate;
        
        % Loop over the number of APs
        for ap = 1:numAPs
            
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
                else % AP to STA
                    % Generate a random secure HE-LTF sequence for the exchange
                    cfgAP{ap}.User{1}.SecureHELTFSequence = dec2hex(randsrc(1,10,(0:15)))';
                    cfg = cfgAP{ap}; % For DL
                end
                
                % Set different channel for UL and DL, assuming that the channel is not reciprocal
                reset(chan)
                
                % Generate HE Ranging NDP transmission
                tx = heRangingWaveformGenerator(cfg);
                
                % Introduce time delay (fractional and integer) in the transmit waveform
                txDelay = heDelaySignal(tx,sampleDelay(ap));
                
                % Pad signal and pass through multipath channel
                txMultipath = chan([txDelay;zeros(50,cfg.NumTransmitAntennas)]);
                
                % Pass waveform through AWGN channel
                rx = awgn(txMultipath,snrVal);
                
                % Perform synchronization and channel estimation
                [chanEstActiveSC,integerOffset] = heRangingSynchronize(rx,cfg);
                
                % Estimate the transmission time between UL and DL
                if ~isempty(chanEstActiveSC) % If packet detection is successful
                    
                    % Estimate fractional delay with MUSIC super-resolution
                    fracDelay = heRangingTOAEstimate(chanEstActiveSC,ofdmInfo.ActiveFFTIndices, ...
                                                     ofdmInfo.FFTLength,sampleRate,numPaths);
                    
                    integerOffset = integerOffset - chDelay; % Account for channel filter delay
                    intDelay = integerOffset/sampleRate; % Estimate integer time delay
                    txTime(l) = intDelay + fracDelay; % Transmission time
                    
                else % If packet detection fails
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
                
                % Estimate the distance between the STA and AP
                distEst(ap,iter,isnr) = (rtt/2)*speedOfLight;
                % Accumulate error to MAE
                rangingError = rangingError + abs(distanceAllAPs(ap) - distEst(ap,iter,isnr));

            else % If packet detection fails
                distEst(ap,iter,isnr) = NaN;
                failedPackets = failedPackets + 1;
            end
        end
    end

    mae = rangingError/((numAPs*numIterations) - failedPackets); % MAE for successful packets
    per(isnr) = failedPackets/(numAPs*numIterations); % PER
    if(per(isnr) > 0.01) % Use only successful packets for ranging and positioning
        warning('wlan:discardPacket','At SNR = %d dB, %d%% of packets were discarded',snrRange(isnr),100*per(isnr));
    end
    disp(['At SNR = ',num2str(snrRange(isnr)),' dB, ','Ranging mean absolute error = ',num2str(mae), ' meters.'])
end

% Reshape to consider all packets within one SNR point as one dataset
rangingError = reshape(abs(distance - distEst),[numAPs*numIterations,numSNR]);
hePlotErrorCDF(rangingError,snrRange) % Cumulative Distribution Function
xlabel('Absolute ranging error (m)');
title('Ranging Error CDF');

%% Plot STA & APs
plotIter = numIterations; plotSNR = numSNR;
plotAP = positionAP(:,:,plotIter,plotSNR); plotSTA = positionSTA(:,plotIter,plotSNR);
figure;
hePlotNodePositions(plotAP, plotSTA); % Plot AP & STA locations

% Plot Trilateration Circles
numAPs = size(positionAP,2);
angles = 0:2*pi/720:2*pi;
for i = 1:numAPs
    x = distEst(i, plotIter, plotSNR) * cos(angles) + plotAP(1,i);
    y = distEst(i, plotIter, plotSNR) * sin(angles) + plotAP(2,i);
    plot(x,y,'k--','LineWidth',1);
    hold on; grid on;
end
title(['Node positions at SNR ' num2str(snrRange(numSNR)) ' dB for iteration #' num2str(plotIter) ])
legend({'AP position','STA poisition','Trilateration circles'},'Location','best','FontSize',10);
legend('boxoff')