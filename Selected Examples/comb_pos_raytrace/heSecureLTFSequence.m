function [HELTF,k,s1,s2,bitsPerSym] = heSecureLTFSequence(cfg,varargin)
%heSecureLTFSequence Secure HE-LTF symbol parameters
%   [HELTF,k,S1,S2,BITSPERSYM] = heSecureLTFSequence(CFG) returns these
%   secure HE-LTF symbol parameters.
% 
%   - HELTF: A matrix of size N-by-NHELTF containing the secure HE-LTF
%     sequence mapped to the corresponding subcarrier indices. N is the
%     number of subcarriers for the corresponding bandwidth and NHELTF is
%     the number of HE-LTF symbols.
%
%   - K: A column vector of corresponding subcarrier indices
%
%   - S1 and S2: A column vector of size 2^P/2-by-NHELTF containing the
%     secure HE-LTF sequence as defined in IEEE P802.11az/D2.0, Equation
%     27-uu and 27-vv, respectively. Number of bits, P are
%     bandwidth-dependent for more information, see IEEE P802.11az/D2.0,
%     Section 27.3.17c.
%
%   - BITSPERSYM: A matrix of size (4P+3)-by-NHELTF of type int8, where
%     NHELTF is the number of HE-LTF symbols of all users in HE-LTF field.
%
%   cfg is a format configuration object of type <a href="matlab:help('heRangingConfig')">heRangingConfig</a>.
%
%   [...] = heSecureLTFSequence(...,UserNumber) returns the secure HE-LTF
%   sequence mapped to the corresponding subcarrier indices, K, S1, S2 and
%   BITSPERSYM for the given UserNumber.

%   Copyright 2020-2023 The MathWorks, Inc.

%#codegen

narginchk(1,2);
validateattributes(cfg,{'heRangingConfig'},{'scalar'},mfilename,'Configuration object');

validateConfig(cfg);
numUsers = numel(cfg.User);
startUserIdx = 1;
endUserIdx = numUsers;

if nargin>1
    startUserIdx = varargin{1};
    if ~isnumeric(startUserIdx)
        error('Second input argument must be numeric scalar.');
    end

    if startUserIdx>numUsers
        error('The UserNumber must not be greater than the total number of users in HE-LTF field.');
    end
end

heRangingInfo = heRangingLTFInfo(cfg);
if nargin>1
    NHELTF = heRangingInfo.NHELTFWithRepetition(startUserIdx);
    heltfSymOffset = zeros(1,numUsers);
    endUserIdx = startUserIdx;
else
    NHELTF = heRangingInfo.NHELTF;
    heltfSymOffset = [0 cumsum(heRangingInfo.NHELTFWithRepetition)];
end

% Reshape the input secure sequence per user into bits per symbols for all users
bitsPerSym = coder.nullcopy(zeros(heRangingInfo.NumSecureHELTFBitsPerSymbol,NHELTF));
octetLength = 8; % 1 octet = 8 bits
for u=startUserIdx:endUserIdx
    % Convert row vector to column of octets
    columnOctets = reshape(cfg.User{u}.SecureHELTFSequence,2,[])';
    % Converting hexadecimal format octets to decimal format
    decOctets = hex2dec(columnOctets);
    bitSeq = reshape(int2bit(decOctets',octetLength,false),[],1); % Convert row vector to column of octets
    repeatFactor = ceil(heRangingInfo.NumSecureHELTFBitsPerUser(u)/numel(bitSeq)); % Sequence repeat factor
    bitSeqExt = [bitSeq; repmat(bitSeq,repeatFactor,1)]; % Extend input sequence
    bitsPerSym(:,heltfSymOffset(u)+(1:heRangingInfo.NHELTFWithRepetition(u))) =  ...
        reshape(bitSeqExt(1:heRangingInfo.NumSecureHELTFBitsPerUser(u)),heRangingInfo.NumSecureHELTFBitsPerSymbol,heRangingInfo.NHELTFWithRepetition(u));
end

% Generate secure HE-LTF sequence as defined in IEEE P802.11az/D2.0,
% Section 27.3.17c.
P = heRangingInfo.P;
[s1,s2] = generateSequence(bitsPerSym,NHELTF,P);
cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
% Map secure HE-LTF symbols into subcarriers
[~,k] = wlan.internal.heLTFSequence(cbw,cfg.HELTFType); % Subcarrier indices for HE-LTF sequence
HELTF = complex(zeros(numel(k),NHELTF));
switch cbw
    case 20
        HELTF(1:2:121,:) = s1(3:63,:); % [-122:2:-2]
        HELTF(125:2:245,:) = s2(2:62,:); % [2:2:122]
    case 40
        HELTF(1:2:241,:) = s1(5:125,:); % [-244:2:-4]
        HELTF(249:2:489,:) = s2(4:124,:); % [4:2:244]
    case 80
        HELTF(1:2:497,:) = s1(5:253,:); % [-500:2:-4]
        HELTF(505:2:1001,:) = s2(4:252,:); % [4:2:500]
    otherwise % 160
        % Add zeros to match CBW160 mapping pattern as defined in IEEE P802.11ax/D4.1, Section 27.3.10.10
        HELTF(1:2:1001,:) = [s1(5:253,:); zeros(3,NHELTF); s1(260:508,:)]; % [-1012:2:-12]
        HELTF(1025:2:2025,:) = [s2(5:253,:); zeros(3,NHELTF); s2(260:508,:)]; % [12:2:1012]
end

% Cyclic shift delay in number of samples. Ts is in ns (1e-9) and chBW
% is in MHz, therefore multiply by 1e-3 to give shift in samples
csInSamples = heRangingInfo.Ts*cbw*1e-3;
% Cylic shift values for all HE-LTF symbols as defined in IEEE
% IEEE P802.11az/D2.0, Equation 27-rr
pBits = bitsPerSym(1:P,1:NHELTF);
cycShift = csInSamples.*sum(pBits.*(2.^(0:P-1).'),1);
% Cyclic shift is applied per HE-LTF symbol
Nfft = 256*cbw/20;
phaseShift = exp(-1i*2*pi*cycShift.*k/Nfft);
HELTF = HELTF.*phaseShift;
bitsPerSym = int8(bitsPerSym);

end

function [s1,s2] = generateSequence(bits,NHELTF,P)
%generateSequence Generate Secure HE-LTF sequence

s1 = coder.nullcopy(complex(zeros(2^P/2,NHELTF)));
s2 = coder.nullcopy(complex(zeros(2^P/2,NHELTF)));
constVal = 2.^[0; 1; 2];
for num=1:NHELTF
    symBits = bits(:,num);
    % Initial value for s1_pm1(0) and s2_pm1(0)
    s1_pm1 = exp(1i*pi/4*(sum(symBits(P:P+2,1).*constVal,1)));
    s2_pm1 = exp(1i*pi/4*(sum(symBits(P+3:P+5,1).*constVal,1)));
    coder.varsize('s1_pm1',[1 512]); % Set to the maximum size 2^P-1 (CBW160)
    coder.varsize('s2_pm1',[1 512]);
    for p=1:P-1
        phi_p = exp(1i*pi/4*(sum(symBits(P+3*p+3:P+3*p+5,1).*constVal)));
        s1_p = [s1_pm1 s2_pm1];
        s2_p = [phi_p*s1_pm1 -phi_p*s2_pm1];
        s1_pm1 = s1_p;
        s2_pm1 = s2_p;
    end
    s1(:,num) = s1_pm1.';
    s2(:,num) = s2_pm1.';
end
end