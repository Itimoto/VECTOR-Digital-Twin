function chanEst = heLTFRangingChannelEstimate(x,cfg,varargin)
%heLTFRangingChannelEstimate Channel estimation using HE-LTF
%   CHANEST = heLTFRangingChannelEstimate(X,CFG) returns the estimated
%   channel between all space-time streams and receive antennas using
%   HE-LTF of an HE ranging NDP. The channel estimates do not include the
%   effect of cyclic shifts at the transmitter.
%
%   CHANEST is an array characterizing the estimated channel for the data
%   and pilot subcarriers. CHANEST is a complex Nst-by-Nsts-by-Nr array
%   characterizing the estimated channel, where Nst is the number of
%   occupied subcarriers, Nsts is the total number of space-time streams,
%   and Nr is the number of receive antennas.
%
%   X is a complex Nst-by-Nsym-by-Nr array containing demodulated
%   concatenated HE-LTF. Nsym is the number of demodulated HE-LTF symbols.
%
%   CFG is a format configuration object of type <a href="matlab:help('heRangingConfig')">heRangingConfig</a>.
%
%   CHANEST = heLTFRangingChannelEstimate(X,CFG,'USERNUMBER') returns the
%   channel estimate for the user of interest. The UserNumber must be a
%   scalar between 1 and 64, inclusive. The default is 1.
%
%   Example:
%   %   Estimate the channel for each user in an HE Ranging NDP
%   %   transmission with a fading channel model.
%
%       % Create packet configuration for two users
%       numUsers = 2;
%       cfg = heRangingConfig(numUsers);
%       cfg.User{1}.SecureHELTFSequence = 'ab';
%       cfg.User{2}.SecureHELTFSequence = 'cd';
%
%       % Generate an HE Ranging NDP waveform
%       tx = heRangingWaveformGenerator(cfg);
%
%       % Get HE-LTF field indices
%       indLTF = heRangingFieldIndices(cfg,'HE-LTF');
%
%       % Process all users
%       for u=1:numUsers
%           % Add channel and noise
%           snr = 20;
%           channel = wlanTGaxChannel;
%           channel.NumTransmitAntennas = cfg.NumTransmitAntennas;
%           channel.NumReceiveAntennas = 1;
%           channel.SampleRate = wlanSampleRate(cfg.ChannelBandwidth);
%           rx = awgn(channel([tx; zeros(10,1)]),snr);
%
%           % Synchronize
%           offset = wlanSymbolTimingEstimate(rx,cfg.ChannelBandwidth);
%           rx = rx(offset:end,:);
%
%           % Extract and OFDM demodulate the HE-LTF for all users
%           rxHETLF = rx(indLTF(u,1): indLTF(u,2),:);
%           demod = heRangingDemodulate(rxHETLF,cfg);
%
%           % Channel estimate for each user
%           chanEst = heLTFRangingChannelEstimate(demod,cfg,u);
%       end
%
%   See also heRangingConfig

%   Copyright 2020-2023 The MathWorks, Inc.

%#codegen

narginchk(2,3);
% Validate the format configuration object is a valid type
validateattributes(cfg,{'ranging.heRangingConfig'},{'scalar'},mfilename,'format configuration object');

userNumber = 1;
if nargin>2
    userNumber = varargin{1};
    validateattributes(userNumber,{'numeric'},{'scalar','>=',1,'<=',64});
end

% Get allocation information
allocInfo = ruInfo(cfg);
if userNumber>numel(cfg.User)
    error('The UserNumber must be less than or equal to the number of users in heRangingConfig object');
end
numSTSUser = cfg.User{userNumber}.NumSpaceTimeStreams;

% Validate symbol type
validateattributes(x,{'double'},{'3d'},mfilename,'HE-LTF OFDM symbol(s)');
[numST,numLTF,numRx] = size(x);

tac = wlan.internal.heRUToneAllocationConstants(allocInfo.RUSizes);
coder.internal.errorIf(numST~=tac.NST,'wlan:wlanChannelEstimate:IncorrectNumSC',tac.NST,numST);
ofdmInfo = wlanHEOFDMInfo('HE-LTF',cfg.ChannelBandwidth,cfg.GuardInterval,[allocInfo.RUSizes allocInfo.RUIndices]);
if numLTF==0
    chanEst = zeros(numST,numSTSUser,numRx);
    return;
end

% Validate the number of OFDM symbols
hezInfo = ranging.heRangingLTFInfo(cfg);
minNumLTF = hezInfo.NHELTFWithRepetition(userNumber);
coder.internal.errorIf(numLTF<minNumLTF,'wlan:he:InvalidNumLTF',numLTF,minNumLTF);
x = x(:,1:minNumLTF,:); % Process the valid demodulated HE-LTF symbols on all receive antennas
N_HE_LTF_Mode = 2;

% Get the HE-LTF sequence
cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
if cfg.SecureHELTF
    % Get HE-LTF secure sequence. IEEE P802.11az/D2.0, Section 27.3.17c
    [HELTF,kHELTFSeq] = ranging.heSecureLTFSequence(cfg,userNumber); % HELTF is of size kHELTFSeq-by-Nltf
else
    % Get HE-LTF sequence. IEEE P802.11ax/D4.1, Section 27.3.10.10
    [HELTF,kHELTFSeq] = wlan.internal.heLTFSequence(cbw,N_HE_LTF_Mode); % HELTF is of size kHELTFSeq-by-1
end

% Extract the RU of interest from the full-bandwidth HELTF
kRU = ofdmInfo.ActiveFrequencyIndices;
[~,ruIdx] = intersect(kHELTFSeq,kRU);
HELTFRU = HELTF(ruIdx,:);

if numSTSUser==1
    % Single STS

    % When more than one LTF we can average over the LTFs for data and
    % pilots to improve the estimate. As there is only one space-time
    % stream, the pilots and data essentially both use the P matrix which
    % does not change per space-time stream (only per symbol), therefore
    % this "MIMO" estimate performs the averaging of the number of symbols.
    chanEst = mimoChannelEstimate(x,HELTFRU,numSTSUser,hezInfo.NHELTFWithoutRepetition(userNumber));

    % Interpolate if HE-LTF compression used
    chanEst = chanEstInterp(chanEst,cbw,N_HE_LTF_Mode,allocInfo.RUSizes,allocInfo.RUIndices);
else
    % MIMO channel estimation as per Perahia, Eldad, and Robert Stacey.
    % Next Generation Wireless LANs: 802.11 n and 802.11 ac. Cambridge
    % University Press, 2013, page 100, Equation 4.39. Remove orthogonal
    % sequence across subcarriers (if used).
    if cfg.SecureHELTF
        % Perform channel estimate on active subcarriers as the same P
        % matrix is used for data and pilot subcarriers.
        chanEst = mimoChannelEstimate(x,HELTFRU,numSTSUser,hezInfo.NHELTFWithoutRepetition(userNumber));
        chanEst = chanEstInterp(chanEst,cbw,N_HE_LTF_Mode,allocInfo.RUSizes,allocInfo.RUIndices);
    else
        % Only perform channel estimate for non-pilot subcarriers as pilots are single stream
        mimoInd = ofdmInfo.DataIndices;
        kMIMO = kRU(mimoInd); % Only data subcarriers MIMO estimates
        chanEst = mimoChannelEstimate(x(mimoInd,:,:),HELTFRU(mimoInd),numSTSUser,hezInfo.NHELTFWithoutRepetition(userNumber));
        % Undo cyclic shift for each STS before averaging and interpolation
        nfft = (cbw/20)*256;
        csh = wlan.internal.getCyclicShiftVal('VHT',numSTSUser,cbw);
        chanEst = wlan.internal.cyclicShiftChannelEstimate(chanEst,-csh,nfft,kMIMO);
        chanEst = chanEstInterp(chanEst,cbw,N_HE_LTF_Mode,allocInfo.RUSizes,allocInfo.RUIndices,mimoInd);
    end
end

end

function chanEstRUInterp = chanEstInterp(chanEstRU,cbw,N_HE_LTF_Mode,ruSize,ruIndex,varargin)
    % Interpolate over pilot locations and compressed subcarriers

    % Get the subcarrier indices within the FFT for the channel estimate
    % input
    Nfft = 256*cbw/20;
    kAct = wlan.internal.heRUSubcarrierIndices(cbw,ruSize,ruIndex)+Nfft/2+1;
    % If the channelEstRU is not the entire RU, then we need to make sure
    % we know the subcarrier indices, so use the ruInd input. For example
    % this allows us to pass in only the data subcarriers.
    if nargin>5
        ruInd = varargin{1};
        kChanEstInputs = kAct(ruInd);
    else
        % Assume chanEstRU is the whole RU
        kChanEstInputs = kAct;
    end

    % Get the indices within the FFT which contain actual estimates
    % (excluding the guard bands). This is how the pattern is structured
    kAll = 1:N_HE_LTF_Mode:Nfft; 

    % Find the subcarrier indices within the FFT which contain actual data
    % within the channel estimate input (kToInterp) and the indices of
    % these within the chanEstDataRU input array (toInterpInd)
    [kToInterp,toInterpInd] = intersect(kChanEstInputs,kAll);

    % Interpolate and extrapolate over all RU subcarrier indices to
    % interpolate over compressed region and pilots
    magPart = abs(chanEstRU(toInterpInd,:,:));
    phasePart = unwrap(angle(chanEstRU(toInterpInd,:,:)));
    combInterp = interp1(kToInterp.',cat(4,magPart,phasePart),kAct);

    % Convert mag and phase to complex
    magPartInterp = combInterp(:,:,:,1);
    phasePartInterp = combInterp(:,:,:,2);
    chanEstRUInterp = magPartInterp.*exp(1i*phasePartInterp);
end

function chanEst = mimoChannelEstimate(x,HELTF,numSTS,numHELTFWithoutRep)
%mimoChannelEstimate MIMO channel estimation
%
%   CHANEST = mimoChannelEstimate(X,HELTF,NUMSTS,NUMHELTFWITHOUTREP)
%   returns channel estimate for each subcarrier given received symbols X,
%   reference sequence HE-LTF, number of space-time streams NUMSTS and
%   number of HE-LTF symbols without repetition.
%
%   X are the received symbols of size Nst-by-Nltf-by-Nrx where NST is the
%   number of subcarriers, Nltf is the number of LTF symbols and Nrx is the
%   number of receive antennas. Note all subcarriers in X must have been
%   modulated with the orthogonal mapping matrix P.
%
%   HELTF is the reference sequence and is sizes Nst-by-1.
%
%   NUMSTS is the number of space-time streams.
%
%   NUMHELTFWITHOUTREP is the number of HE-LTF symbols in HE-LTF excluding
%   repetition of HE-LTF symbols.
%
%   MIMO channel estimation as per Perahia, Eldad, and Robert Stacey.
%   Next Generation Wireless LANs: 802.11 n and 802.11 ac. Cambridge
%   university press, 2013, page 100, Eq 4.39.

[Nst,Nltf,Nrx] = size(x);
numRep = Nltf/numHELTFWithoutRep; % Number of HE-LTF repetitions
% Get P matrix
P = wlan.internal.mappingMatrix(numHELTFWithoutRep);
Puse = P(1:numSTS,1:numHELTFWithoutRep)'; 
PRep = repmat(Puse,numRep,1);
denom = Nltf.*HELTF;
chanEst = coder.nullcopy(zeros(Nst,numSTS,Nrx,'like',x));
for k = 1:Nrx
    rxsym = squeeze(x(:,(1:Nltf),k)); % Symbols on 1 receive antenna
    for l = 1:numSTS
        chanEst(:,l,k) = rxsym./denom*PRep(:,l);
    end
end

end
