function demod = heRangingDemodulate(rx,cfg,varargin)
%heRangingDemodulate Demodulate HE-LTF
%   DEMOD = heRangingDemodulate(RX,CFG) demodulates the time-domain HE-LTF
%   signal RX using OFDM demodulation parameters.
%
%   DEMOD is the demodulated frequency-domain signal, returned as a complex
%   matrix or 3-D array of size Nst-by-Nsym-by-Nr. Nst is the number of
%   active (occupied) subcarriers in the field. Nsym is the number of OFDM
%   symbols. Nr is the number of receive antennas.
%
%   RX is the received time-domain signal, specified as a complex matrix of
%   size Ns-by-Nr, where Ns represents the number of time-domain samples.
%   If Ns is not an integer multiple of the OFDM symbol length for the
%   specified field, then mod(Ns,symbol length) trailing samples are
%   ignored.
%
%   CFG is a format configuration object of type <a href="matlab:help('heRangingConfig')">heRangingConfig</a>.
%
%   DEMOD = heRangingDemodulate(...,'OFDMSymbolOffset',SYMOFFSET) specifies
%   the optional OFDM symbol sampling offset as a fraction of the guard
%   interval length between 0 and 1, inclusive. When unspecified, a value
%   of 0.75 is used.
%
%   Examples:
%
%   Example 1:
%   % Demodulate the HE-LTF field of an HE Ranging waveform for a single
%   % user.
%
%      % Generate waveform
%      cfg = heRangingConfig;
%      tx = heRangingWaveformGenerator(cfg,'WindowTransitionTime',0);
%
%      % Generate field indices and extract HE-LTF
%      ind = heRangingFieldIndices(cfg,'HE-LTF');
%      heltf = tx(ind(1):ind(2),:);
%
%      % Demodulate HE-LTF field
%      demodLTF = heRangingDemodulate(heltf,cfg);  
%
%   Example 2:
%   % Demodulate the HE-LTF field of an HE Ranging waveform for two users.
%
%      % Configure users
%      cfg = heRangingConfig(2);
%      cfg.User{1}.NumHELTFRepetition = 2;
%      cfg.User{1}.SecureHELTFSequence = 'a2b3';
%      cfg.User{2}.NumHELTFRepetition = 3;
%      cfg.User{2}.SecureHELTFSequence = 'c4d5';
%
%      % Generate waveform
%      tx = heRangingWaveformGenerator(cfg,'WindowTransitionTime',0);
%
%      % Generate field indices and extract HE-LTF
%      ind = heRangingFieldIndices(cfg,'HE-LTF');
%
%      % Demodulate HE-LTF field for user-1
%      heltf = tx(ind(1,1):ind(1,2),:);
%      demodLTFUser1 = heRangingDemodulate(heltf,cfg); 
%
%      % Demodulate HE-LTF field for user-2
%      heltf = tx(ind(2,1):ind(2,2),:);
%      demodLTFUser2 = heRangingDemodulate(heltf,cfg); 
%
%   See also heLTFRangingChannelEstimate

%   Copyright 2020-2022 The MathWorks, Inc.

%#codegen

narginchk(2,4);

% Validate the format configuration object is a valid type
validateattributes(cfg,{'heRangingConfig'},{'scalar'},mfilename,'format configuration object');

% Validate rx and field name
validateattributes(rx,{'double'},{'2d','finite','nonempty'},mfilename,'rx');
[numSamples,numRx] = size(rx);

infoRU = ruInfo(cfg);
cfgOFDM = wlanHEOFDMInfo('HE-LTF',cfg.ChannelBandwidth,cfg.GuardInterval,[infoRU.RUSizes infoRU.RUIndices]);

nvp = wlan.internal.demodNVPairParse(varargin{:});
uncompFFTLength = cfgOFDM.FFTLength;
N_HE_LTF_Mode = cfg.HELTFType;
compFFTLength = uncompFFTLength/N_HE_LTF_Mode;

% Calculate number of symbols to demodulate
Ns = compFFTLength+cfgOFDM.CPLength;

% Validate input length
validateMinInputLength(numSamples,Ns,cfgOFDM.CPLength,cfg.SecureHELTF);

if numSamples==0
    demod = complex(zeros(cfgOFDM.NumTones,0,numRx)); % Return empty for 0 samples
    return;
end
    numSym = floor(numSamples/Ns);
    numSamples = numSym*(compFFTLength+cfgOFDM.CPLength);
    % OFDM demodulate
    prmStr = struct;
    prmStr.NumReceiveAntennas = numRx;
    prmStr.FFTLength = compFFTLength; % Use compressed FFT length
    prmStr.NumSymbols = numSym;
    prmStr.SymbolOffset = round(nvp.SymOffset*cfgOFDM.CPLength);
    prmStr.CyclicPrefixLength = cfgOFDM.CPLength;
    compressedfftout = coder.nullcopy(complex(zeros(compFFTLength,numSym,numRx))); % For codegen
    if cfg.SecureHELTF
        numSamples = numSamples+cfgOFDM.CPLength; % Account for trailing guard interval after last OFDM symbol
        compressedfftout(:) = zpDemodulate(rx(1:numSamples,:),prmStr);
    else
        compressedfftout(:) = comm.internal.ofdm.demodulate(rx(1:numSamples,:),prmStr);
    end

    % Decompress HE-LTF
    fftout = complex(zeros(uncompFFTLength,numSym,numRx));
    fftout(1:N_HE_LTF_Mode:end,:,:) = compressedfftout;

    % Extract active subcarriers from full FFT
    demod = fftout(cfgOFDM.ActiveFFTIndices,:,:);

    % Scale by number of tones and FFT length at transmitter
    demod = demod*sqrt(cfgOFDM.NumTones/N_HE_LTF_Mode)/compFFTLength;
end

function validateMinInputLength(numSamples,Ns,cpLen,secureHELTF)
% Validate number of input samples

    if secureHELTF
        Ns = Ns+cpLen;
    end

    minInputLength = sum(Ns);
    coder.internal.errorIf(numSamples>0 && numSamples<minInputLength,'wlan:shared:ShortDataInput',minInputLength);
end

function y = zpDemodulate(x,prmStr)
%   Y = zpDemodulate(X,PRMSTR) performs OFDM demodulation on input X, using
%   the parameters specified in the structure PRMSTR.
%
%   Y is the demodulated frequency-domain signal, returned as a complex
%   matrix or 3-D array of size Nst-by-Nsym-by-Nr. Nst is the number of
%   subcarriers in the field. Nsym is the number of OFDM symbols. Nr is the
%   number of receive antennas.
%
%   X is the time-domain signal, specified as a complex matrix of size
%   Ns-by-Nr, where Ns represents the number of time-domain samples.
%
%   PRMSTR is a structure with the following fields:
%       FFTLength
%       CyclicPrefixLength
%       NumSymbols
%       SymbolOffset
%       NumReceiveAntennas

    fftLen = prmStr.FFTLength;
    giLen  = prmStr.CyclicPrefixLength(1);
    symOffset = prmStr.SymbolOffset;
    symLength = fftLen+giLen; % OFDM Symbol length in samples
    NSYM = prmStr.NumSymbols;
    NRX = size(x,2);

    % Discard intial guard interval accounting for OFDM symbol offset
    todemod = reshape(x(symOffset+1:(end-(giLen-symOffset)),:),[symLength,NSYM,NRX]);

    % Overlap start of the symbol with the next guard interval
    todemod(1:giLen,:,:) = todemod(1:giLen,:,:)+todemod(fftLen+(1:giLen),:,:);

    % Remove GI
    todemod = todemod(1:fftLen,:,:);

    % Move symbol offset samples to end of symbol to maintain phase
    todemod = [todemod((giLen-symOffset)+1:fftLen,:,:); todemod(1:(giLen-symOffset),:,:)];

    % FFT
    postFFT = fft(todemod,[],1);
    y = fftshift(postFFT,1); 
end