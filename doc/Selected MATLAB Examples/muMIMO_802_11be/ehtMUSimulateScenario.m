function throughput = ehtMUSimulateScenario(cfg,tgax,cfgSim)
%ehtMUSimulateScenario Simulate the EHT MU scenario
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   THROUGHPUT = ehtMUSimulateScenario(CFG,TGAX,CFGSIM) simulates a EHT MU
%   scenario and returns the throughput.
%
%   CFG is the configuration of the EHT MU transmission and is a format
%   configuration object of type
%   <a href="matlab:help('wlanEHTMUConfig')">wlanEHTMUConfig</a>.
%
%   TGAX is a cell array containing channels models between the transmitter
%   and each user.
%
%   CFGSIM is a structure containing the simulation parameters.

%   Copyright 2022 The MathWorks, Inc.

% Simulation options
numPackets = cfgSim.NumPackets;
pathloss = cfgSim.Pathloss;
noiseFloor = cfgSim.NoiseFloor;
idleTime = cfgSim.IdleTime; % us
transmitPower = cfgSim.TransmitPower; % dBm

% Perform synchronization with 11ac components
timingSync = true;
cfoCorrection = true;
pilotPhaseTracking = true;

% Get the field indices for extract fields from the PPDU
ind = wlanFieldIndices(cfg);

% Get allocation info
allocInfo = ruInfo(cfg);

s = RandStream.getGlobalStream(); % Save random stream to restore at end
numRUs = allocInfo.NumRUs;
numPL = numel(pathloss);
packetErrorRate = zeros(numPL,allocInfo.NumUsers);
throughput = zeros(numPL,allocInfo.NumUsers);
numPacketErrors = zeros(numPL,allocInfo.NumUsers);
fs = wlanSampleRate(cfg);
A = 10.^((transmitPower-30)/20); % Voltage gain (attenuation)
chanBW = cfg.ChannelBandwidth;

% Simulate all path losses
for ipl = 1:numPL

    % Clone the channel before running to ensure the same realization is
    % used for all path losses and schemes
    tgaxClone = cloneChannels(tgax);

    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = ipl;
    RandStream.setGlobalStream(stream);

    % Create an instance of the AWGN channel per SNR point simulated
    awgnChannel = comm.AWGNChannel('NoiseMethod','Variance',...
        'Variance',10^((noiseFloor-30+pathloss(ipl))/10));

    % Simulate multiple packets
    for pktIdx = 1:numPackets
        % Generate random data for all users
        psdu = psduLength(cfg);
        txPSDU = cell(allocInfo.NumUsers,1);
        for userIdx = 1:allocInfo.NumUsers
            txPSDU{userIdx} = randi([0 1],psdu(userIdx)*8,1,'int8');
        end

        % Generate waveform with idle period
        tx = wlanWaveformGenerator(txPSDU,cfg,'IdleTime',idleTime*1e-6);

        % Scale to achieve desired transmit power
        tx = tx*A;

        % Per user processing
        userIdx = 1;

        % Loop over RUs
        for ruIdx = 1:numRUs
            % Loop over users within an RU
            for ruUserIdx = 1:allocInfo.NumUsersPerRU(ruIdx)
                % Pass waveform through a channel for each user
                rx = tgaxClone{userIdx}(tx);

                % Add noise at receiver (assuming pathloss and noise floor)
                rx = awgnChannel(rx);

                if timingSync
                    % Packet detect and determine coarse packet offset
                    coarsePktOffset = wlanPacketDetect(rx,chanBW);
                    if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
                        numPacketErrors(ipl,userIdx) = numPacketErrors(ipl,userIdx)+1;
                        userIdx = userIdx+1;
                        continue; % Go to next loop iteration
                    end
                else
                   coarsePktOffset = 0;  %#ok<UNRCH>
                end

                if cfoCorrection
                    % Extract L-STF and perform coarse frequency offset correction
                    lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
                    coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
                    rx = frequencyOffset(rx,fs,-coarseFreqOff);
                end

                if timingSync
                    % Extract the non-HT fields and determine fine packet offset
                    nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
                    finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

                    % Determine final packet offset
                    pktOffset = coarsePktOffset+finePktOffset;
                else
                    pktOffset = 4; %#ok<UNRCH>
                end

                % If packet detected outwith the range of expected delays from
                % the channel modeling; packet error
                if pktOffset>50
                    numPacketErrors(ipl,userIdx) = numPacketErrors(ipl,userIdx)+1;
                    userIdx = userIdx+1;
                    continue; % Go to next loop iteration
                end

                if cfoCorrection
                    % Extract L-LTF and perform fine frequency offset correction
                    rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
                    fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
                    rx = frequencyOffset(rx,fs,-fineFreqOff);
                end

                % EHT-LTF demodulation and channel estimation
                rxEHTLTF = rx(pktOffset+(ind.EHTLTF(1):ind.EHTLTF(2)),:);
                ehtltfDemod = wlanEHTDemodulate(rxEHTLTF,'EHT-LTF',cfg,ruIdx);
                [chanEst,pilotEst] = wlanEHTLTFChannelEstimate(ehtltfDemod,cfg,ruIdx);

                % EHT data demodulate
                rxData = rx(pktOffset+(ind.EHTData(1):ind.EHTData(2)),:);
                demodSym = wlanEHTDemodulate(rxData,'EHT-Data',cfg,ruIdx);

                % Pilot phase tracking
                if pilotPhaseTracking
                    % Average single-stream pilot estimates over symbols (2nd dimension)
                    pilotEstTrack = mean(pilotEst,2);
                    demodSym = wlanEHTTrackPilotError(demodSym,pilotEstTrack,cfg,'EHT-Data',ruIdx);
                end

                % Estimate noise power in EHT fields
                ruOFDMInfo = wlanEHTOFDMInfo('EHT-Data',cfg,ruIdx);
                nVarEst = ehtNoiseEstimate(demodSym(ruOFDMInfo.PilotIndices,:,:),pilotEst,cfg,ruIdx);

                % Equalize
                [eqSym,csi] = ehtEqualizeCombine(demodSym,chanEst,nVarEst,cfg,userIdx);

                % Discard pilot subcarriers
                rxDataUser = eqSym(ruOFDMInfo.DataIndices,:,:);
                csiData = csi(ruOFDMInfo.DataIndices,:);

                % Demap and decode bits
                rxPSDU = wlanEHTDataBitRecover(rxDataUser,nVarEst,csiData,cfg,userIdx);

                % Compare bit error
                packetError = ~isequal(txPSDU{userIdx},rxPSDU);
                numPacketErrors(ipl,userIdx) = numPacketErrors(ipl,userIdx)+packetError;
                userIdx = userIdx+1;
            end
        end
    end
    % Calculate packet error rate (PER) at SNR point
    packetErrorRate(ipl,:) = numPacketErrors(ipl,:)/numPackets;

    % Calculate the raw throughput
    txtime = transmitTime(cfg,'microseconds'); % Packet duration
    for userIdx = 1:allocInfo.NumUsers
        throughput(ipl,userIdx) = 1e-6*((numPackets-numPacketErrors(ipl,userIdx))*(8*cfg.User{userIdx}.APEPLength)/(numPackets*(txtime+idleTime)*1e-6));
    end

    disp([' Pathloss ' num2str(pathloss(ipl),'%2.1f') ' dB,'...
        ' AP throughput ' num2str(sum(throughput(ipl,:)),'%2.1f') ' Mbps']);
end

RandStream.setGlobalStream(s); % Restore random stream
end

function tgaxClone = cloneChannels(tgax)
% Clone the channels before running to ensure the same realization can be
% used in the next simulation
    tgaxClone = cell(size(tgax));
    numUsers = numel(tgax);
    for i=1:numUsers
        tgaxClone{i} = clone(tgax{i});
    end
end


