function y = heRangingLTF(cfg)
%heRangingLTF HE Ranging NDP Long Training Field (HE-LTF)
%   Y = heRangingLTF(CFG) generates an HE Ranging Long Training Field
%   (HE-LTF) time-domain signal for the HE Ranging NDP format as defined in
%   IEEE P802.11az/D2.0 standard.
%
%   Y is the time-domain HE-LTF signal. It is a complex matrix of size
%   Ns-by-Nt where Ns represents the number of time-domain samples and Nt
%   represents the number of transmit antennas.
%
%   CFG is the format configuration object of type <a href="matlab:help('heRangingConfig')">heRangingConfig</a>
%
%   % Example: Generate an HE Ranging LTF field for a 80MHz.
%
%   cfg = heRangingConfig('ChannelBandwidth','CBW80');
%   y = heRangingLTF(cfg);
%   plot(abs(y));

%   Copyright 2020-2023 The MathWorks, Inc.

%#codegen

validateattributes(cfg,{'ranging.heRangingConfig'},{'scalar'},mfilename,'Configuration object');
validateConfig(cfg);
numTx = cfg.NumTransmitAntennas;
cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
numUsers = numel(cfg.User);      % Number of user in HE-LTF field
allocationInfo = ruInfo(cfg);    % RU allocation info
hezInfo = ranging.heRangingLTFInfo(cfg); % Get HE-LTF field information
Nltf = hezInfo.NHELTF;           % Number of HE-LTF symbols
maxNumSpaceTimeStream = cfg.NumSpaceTimeStreams; % Max number of space-time streams between all users

if cfg.SecureHELTF
    % Get HE-LTF secure sequence. IEEE P802.11az/D2.0, Section 27.3.17c
    [HELTF,kHELTFSeq] = ranging.heSecureLTFSequence(cfg); % HELTF is of size kHELTFSeq-by-Nltf
else
    % Get HE-LTF sequence. IEEE P802.11ax/D4.1, Section 27.3.10.10
    [HELTF,kHELTFSeq] = wlan.internal.heLTFSequence(cbw,cfg.HELTFType); % HELTF is of size kHELTFSeq-by-1
end

Nfft = 256*cbw/20;
ofdmGrid = complex(zeros(Nfft,Nltf,numTx));
kRUPuncture = wlan.internal.hePuncturedRUSubcarrierIndices(cfg);

% Extract the RU of interest from the full-bandwidth HELTF
kRUFull = wlan.internal.heRUSubcarrierIndices(cbw,allocationInfo.RUSizes(1),allocationInfo.RUIndices(1));
if ~isempty(kRUPuncture)
    kRU = setdiff(kRUFull,kRUPuncture); % Discard punctured subcarriers
else
    kRU = kRUFull;
end
[~,ruIdx] = intersect(kHELTFSeq,kRU);
HELTFRU = HELTF(ruIdx,:);

% Create mapping matrix for all HE-LTF symbols per user
% Generate P and R matrix for all HE-LTF symbols within an HE-LTF field
PRep = complex(zeros(maxNumSpaceTimeStream,Nltf));
RRep = complex(zeros(maxNumSpaceTimeStream,Nltf));
heltfSymOffset = [0 cumsum(hezInfo.NHELTFWithRepetition)];
numLTF = hezInfo.NHELTFWithoutRepetition;

powerBoostFactor = 1; % No additional power boost (alpha)
cardKr = numel(kRU);
cardKHELTFr = cardKr/cfg.HELTFType;
userSTSScalingFactor = zeros(1,Nltf);
scalingConstant = powerBoostFactor*sqrt(cardKr);

for u = 1:numUsers
    % Mapping matrix for all symbols and users. IEEE P802.11ax/D4.1, Equation 27-52
    numSTS = cfg.User{u}.NumSpaceTimeStreams;
    numRep = cfg.User{u}.NumHELTFRepetition;
    % Calculate P and R values in a matrix
    Pheltf = wlan.internal.mappingMatrix(numLTF(u));
    % Get P matrix for the number of HE-LTF symbols
    P = Pheltf((1:numSTS).',1:numLTF(u));
    % Same P matrix is used for all the repeated HE-LTF symbols
    P = repmat(P,1,numRep);
    PRep(1:numSTS,(1:hezInfo.NHELTFWithRepetition(u))+heltfSymOffset(u)) = P;
    % Get R matrix for the number of space-time streams
    R = repmat(Pheltf(1,1:numLTF(u)),numSTS,1);
    % Same R matrix is used for all the repeated HE-LTF symbols
    R = repmat(R,1,numRep);
    RRep(1:numSTS,(1:hezInfo.NHELTFWithRepetition(u))+heltfSymOffset(u)) = R;

    % STS scaling per user
    % Calculate per RU scaling as per IEEE P802.11ax/D4.1, Section 27.3.9,
    % Equation 27-5. Note we scale by sqrt(cardK)/sqrt(cardKHELTFr). This
    % causes us to normalize the power in the RU. The HE-LTF symbols for
    % the user are scalled by the number of space-time streams for the user
    % of interest.
    ruScalingFactor = scalingConstant/sqrt(cfg.User{u}.NumSpaceTimeStreams*cardKHELTFr);
    userSTSScalingFactor(1,(1:hezInfo.NHELTFWithRepetition(u))+heltfSymOffset(u)) = repmat(ruScalingFactor,1,hezInfo.NHELTFWithRepetition(u));
end

ruGrid = complex(zeros(numel(kRU),Nltf,maxNumSpaceTimeStream));
% Indices of data and pilot subcarriers within the occupied RU
kPilot = wlan.internal.hePilotSubcarrierIndices(cbw,allocationInfo.RUSizes);
ruInd = wlan.internal.heOccupiedSubcarrierIndices(kRU,kPilot);

% Indices of data and pilot subcarriers within the occupied RU
for k = 1:Nltf
    if cfg.SecureHELTF % IEEE P802.11ax/D4.1, Equation 27-52
        % Data and pilots are both multiplied by P matrix so whole RU
        % multiplied by P matrix
        ruGrid(:,k,:) = squeeze(HELTFRU(:,k,:)).*PRep(:, k).';
    else
        % Single stream pilots
        ruGrid(ruInd.Data,k,:) = HELTFRU(ruInd.Data).*PRep(:, k).';
        ruGrid(ruInd.Pilot,k,:) = HELTFRU(ruInd.Pilot).*RRep(:, k).';
    end
end

if ~cfg.SecureHELTF
    % Cyclic shift is applied per space time stream
    ruGrid = wlan.internal.heCyclicShift(ruGrid,cbw,kRU,(1:maxNumSpaceTimeStream).');
end

% Map subcarriers to full FFT grid and scale
n_Scale = 1; % Scaling factor for an HE SU format. IEEE P802.11ax/D4.1, Section 27.3.10.10
ofdmGrid(kRU+Nfft/2+1,:,1:maxNumSpaceTimeStream) = ruGrid.*userSTSScalingFactor*n_Scale;

% Overall scaling factor
allScalingFactor = Nfft/sqrt(sum(powerBoostFactor.^2.*cardKr));

switch cfg.GuardInterval
    case 0.8
        CPLenSamples = 0.8*cbw;
    otherwise % 1.6
        CPLenSamples = 1.6*cbw;
end

if cfg.SecureHELTF % OFDM modulate with zero power GI
    postShift = ifftshift(ofdmGrid(1:2:end,:,:),1);
    postIFFT = ifft(postShift,[],1);
    postCP = [complex(zeros(CPLenSamples,Nltf,numTx)); postIFFT]; % Add zero power CP
    y = reshape(postCP,[(Nfft/2+CPLenSamples)*Nltf numTx]).*allScalingFactor/2; %/2 to account for NFFT/2 actual FFT size
else
    y = wlan.internal.wlanOFDMModulate(ofdmGrid(1:2:end,:,:),CPLenSamples).*allScalingFactor/2; %/2 to account for NFFT/2 actual FFT size
end

end