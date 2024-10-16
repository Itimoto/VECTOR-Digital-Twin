function hezInfo = heRangingLTFInfo(cfg)
%heRangingLTFInfo HE-LTF information of Ranging NDP packet
%   hezInfo = heRangingLTFInfo(cfg) returns a structure containing
%   HE-LTF field information for the format configuration object cfg.
%   The structure hezInfo has the following fields:
%
%   NHELTF                      - Number of HE-LTF symbols in HE-LTF
%
%   NHELTFWithRepetition        - Number of HE-LTF symbols in HE-LTF
%                                 including repetition of HE-LTF
%                                 symbols for all users as a vector of size
%                                 1-by-NumUsers.
%
%   NHELTFWithoutRepetition     - Number of HE-LTF symbols in HE-LTF
%                                 excluding repetition of HE-LTF symbols as
%                                 a vector of size 1-by-NumUsers.
%
%   NumSecureHELTFBitsPerSymbol - Number of secure bits required to encode
%                                 an HE-LTF symbol.
%
%   NumSecureHELTFBitsPerUser   - Number of secure bits required to encode
%                                 the HE-LTF symbols for all users as a
%                                 vector of size 1-by-NumUsers.
%
%   Ts                          - Ts is the cyclic shift in nanosecond as
%                                 defined in IEEE P802.11az/D2.0, Section
%                                 27.3.17b.
%
%   P                           - Number of CSD bits to extract from the
%                                 randomized bits for the given bandwidth
%                                 as defined in IEEE P802.11az/D2.0,
%                                 Section 27.3.17b.
%
%   cfg is a format configuration object of type <a href="matlab:help('heRangingConfig')">heRangingConfig</a>.

%   Copyright 2020 The MathWorks, Inc.

%#codegen

validateattributes(cfg,{'heRangingConfig'},{'scalar'},mfilename,'Configuration object');
cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
switch cbw
    case 20
        P = 7;   % Number of P bits, IEEE P802.11az/D2.0, Section 27.3.17c
        Ts = 50; % Cyclic shift in nsec, IEEE P802.11az/D2.0, Section 27.3.17c
    case 40
        P = 8;
        Ts = 25;
    case 80
        P = 9;
        Ts = 12.5;
    otherwise % 160MHz
        P = 10;
        Ts = 6.25;
end

numUsers = numel(cfg.User);
numHELTFSymPerUserWithRep = coder.nullcopy(zeros(1,numUsers));
numHELTFSymPerUserWithoutRep = coder.nullcopy(zeros(1,numUsers));
for u=1:numUsers
   numHELTFSymPerUserWithoutRep(u) = wlan.internal.numVHTLTFSymbols(cfg.User{u}.NumSpaceTimeStreams);
   numHELTFSymPerUserWithRep(u) = numHELTFSymPerUserWithoutRep(u)*cfg.User{u}.NumHELTFRepetition;
end
numHELTFTotal = sum(numHELTFSymPerUserWithRep);
numBitsPerSym = 4*P+3; % Number of bits required to encode HE-LTF symbol
numBitsPerUser = numBitsPerSym.*numHELTFSymPerUserWithRep; % Total number of bits required to encode HE-LTF user
assert(numBitsPerSym*numHELTFTotal==sum(numBitsPerUser));
hezInfo = struct(...
            'NHELTF',                      numHELTFTotal, ...
            'NHELTFWithRepetition',        numHELTFSymPerUserWithRep, ...
            'NHELTFWithoutRepetition',     numHELTFSymPerUserWithoutRep, ...
            'NumSecureHELTFBitsPerSymbol', numBitsPerSym, ...
            'NumSecureHELTFBitsPerUser',   numBitsPerUser, ...
            'Ts',                          Ts, ... % In nanaseconds
            'P',                           P);
end