function steeringMat = ehtUserBeamformingFeedback(rx,cfgNDP,varargin)
%ehtUserBeamformingFeedback EHT user beamforming feedback
%
%   STEERINGMAT = ehtUserBeamformingFeedback(RX,CFGNDP) returns the
%   steering matrix recommended to beamform towards the user of interest.
%   The steering matrix is calculated using SVD.
%
%   STEERINGMAT = ehtUserBeamformingFeedback(RX,CFGNDP,CFOCOMP) performs
%   carrier frequency offset estimation and compensation on the received
%   waveform before calculating the steering matrix if CFOCOMP is set to
%   TRUE.
%
%   STEERINGMAT is a Nst-by-Nsts-by-Ntx array containing the recommended
%   full band beamforming steering matrix. Nst is the number of occupied
%   subcarriers, Nsts is the maximum number of space-time streams, and Ntx
%   is the number of transmit antennas.
%
%   RX is the received NDP at the station.
%
%   CFGNDP is the format configuration object of type <a href="matlab:help(wlanEHTMUConfig')">wlanEHTMUConfig</a>.

%   Copyright 2022 The MathWorks, Inc.

narginchk(2,3);

if nargin == 3
    cfoCompensate = varargin{1};
    validateattributes(cfoCompensate,{'logical'},{'nonnan','finite'},mfilename,'',3);
else
    cfoCompensate = false;
end

chanBW = cfgNDP.ChannelBandwidth;
ind = wlanFieldIndices(cfgNDP);
fs = wlanSampleRate(cfgNDP);

% Packet detect and determine coarse packet offset
coarsePktOffset = wlanPacketDetect(rx,chanBW);
if isempty(coarsePktOffset) % If empty no L-STF detected; try 0
    % Synchronization failed, return empty steering matrix
    steeringMat = [];
    return;
end

if cfoCompensate
    % Extract L-STF and perform coarse frequency offset correction
    lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
    rx = frequencyOffset(rx,fs,-coarseFreqOff);
end

% Extract the non-HT fields and determine fine packet offset
nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

% Determine final packet offset
pktOffset = coarsePktOffset+finePktOffset;

% If packet detected outwith the range of expected delays from
% the channel modeling; packet error
if pktOffset>50
    % Synchronization failed, return empty steering matrix
    steeringMat = [];
    return;
end

if cfoCompensate
    % Extract L-LTF and perform fine frequency offset correction
    rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
    fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
    rx = frequencyOffset(rx,fs,-fineFreqOff);
end

% Demodulate EHT-LTF with info for all RUs
rxEHTLTF = rx(pktOffset+(ind.EHTLTF(1):ind.EHTLTF(2)),:);
demodEHTLTF = wlanEHTDemodulate(rxEHTLTF,'EHT-LTF',cfgNDP);

% Extract the demodulated EHT-LTF RU
demodEHTLTFRU = demodEHTLTF;

% Channel estimate with info from current RU
chanEstUser = wlanEHTLTFChannelEstimate(demodEHTLTFRU,cfgNDP);

% Get cyclic shift and inverse
numTx = cfgNDP.User{1}.NumSpaceTimeStreams; % The sounding packet is configured for one space-time stream per antenna
csh = -wlan.internal.getCyclicShiftVal('VHT',numTx,wlan.internal.cbwStr2Num(chanBW));
% Indices of active subcarriers in the NDP
ndpOFDMInfo = wlanEHTOFDMInfo('EHT-Data',cfgNDP);
chanEstMinusCSD = permute( ...
    wlan.internal.cyclicShift(permute(chanEstUser,[1 3 2]),csh,ndpOFDMInfo.FFTLength,ndpOFDMInfo.ActiveFrequencyIndices), ...
    [2 3 1]); % Nr-by-Nsts-by-Nst

% Compute the feedback matrix using singular value decomposition for the
% streams allocated to the user
[~,~,V] = pagesvd(chanEstMinusCSD,'econ');
steeringMat = permute(V,[3 2 1]); % Permute to Nst-by-Nsts-by-Ntx

end
