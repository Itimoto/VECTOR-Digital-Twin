function indices = heRangingFieldIndices(cfgFormat,varargin)
%heRangingFieldIndices Generate field indices for an HE Ranging packet
%   INDICES = heRangingFieldIndices(CFGFORMAT) returns the start and end
%   time-domain sample indices for all fields in a packet relative to the
%   first sample.
%
%   INDICES is a structure array of field names for the specified
%   configuration and contains the start and end indices of all fields in a
%   packet.
%
%   CFGFORMAT is a format configuration object of type <a href="matlab:help('heRangingConfig')">heRangingConfig</a>.
%
%   INDICES = heRangingFieldIndices(CFGFORMAT,FIELDNAME) returns the start
%   and end time-domain sample indices for the specified FIELDNAME in a
%   packet. INDICES is a row vector of size two representing the start and
%   end sample indices of the specified field.
%
%   FIELDNAME is a character vector or string specifying the field of
%   interest.
%
%   The field of interest must be 'L-STF', 'L-LTF', 'L-SIG', 'RL-SIG',
%   'HE-SIG-A', HE-SIG-B', 'HE-STF','HE-LTF','HE-Data' or 'HE-PE'.
%
%   When you specify a secure HE-LTF, the function returns the field
%   indices as a matrix of size R-by-2, where R is the number of users.
%   Each row of a matrix represents the start and end indices of the HE-LTF
%   for a user including the zero-power GI after the last HE-LTF symbol.
%
%   Example 1:
%   %   Extract the HE-LTF field for a received HE Ranging transmission.
%
%      cfg = heRangingConfig; % Create packet configuration
%
%      % Generate transmit waveform
%      txWaveform = heRangingWaveformGenerator(cfg);
%
%      % Recover HE-LTF field
%      ind = heRangingFieldIndices(cfg);
%      rxHELTF = txWaveform(ind.HELTF(1):ind.HELTF(2));
%
%   Example 2:
%   %   Extract the HE-LTF field from an HE Ranging packet for multiple
%   %   users.
%
%      % Create packet configuration
%      cfg = heRangingConfig(2);
%
%      % Generate transmit waveform
%      txWaveform = heRangingWaveformGenerator(cfg);
%
%      % Recover HE-LTF for user-2
%      heLTFIdx = heRangingFieldIndices(cfg,'HE-LTF');
%      rxHELTF = txWaveform(heLTFIdx(2,1):heLTFIdx(2,2));

%   Copyright 2020-2021 The MathWorks, Inc.

%#codegen

% Validate input is a class object
validateattributes(cfgFormat,{'heRangingConfig'},{'scalar'},mfilename,'format configuration object');
narginchk(1,2)
if nargin == 1
    indices = getHEIndices(cfgFormat);
else
    fieldType = varargin{1};
    coder.internal.errorIf(~(ischar(fieldType) || (isstring(fieldType) && isscalar(fieldType))) || ... 
            ~any(strcmpi(fieldType,{'L-STF','L-LTF','L-SIG','RL-SIG','HE-SIG-A','HE-SIG-B','HE-STF','HE-LTF','HE-Data','HE-PE'})), ...
            'wlan:wlanFieldIndices:InvalidFieldTypeHE');
    indices = getHEIndices(cfgFormat,fieldType);
end
end

function ind = getHEIndices(cfg,varargin)

    validateConfig(cfg);
    trc = wlan.internal.heTimingRelatedConstants(cfg.GuardInterval,cfg.HELTFType,4,0,0); % Use default preFEC padding factor
    cbw = wlan.internal.cbwStr2Num(cfg.ChannelBandwidth);
    sf = cbw*1e-3; % Scaling factor to convert bandwidth and time in ns to samples

    % L-STF
    nFieldSamples = trc.TLSTF*sf;
    indLSTF = uint32([1 nFieldSamples]);
    numCumSamples = nFieldSamples;

    if nargin>1 && strcmp(varargin{1},'L-STF')
        ind = indLSTF;
        return;
    end

    % L-LTF
    nFieldSamples = trc.TLLTF*sf;
    indLLTF = uint32([numCumSamples+1 numCumSamples+nFieldSamples]);
    numCumSamples = numCumSamples+nFieldSamples;
    nFieldSamples = trc.TLSIG*sf;

    if nargin>1 && strcmp(varargin{1},'L-LTF')
        ind = indLLTF;
        return;
    end

    % L-SIG
    indLSIG = uint32([numCumSamples+1 numCumSamples+nFieldSamples]);
    numCumSamples = numCumSamples+nFieldSamples;
    if nargin>1 && strcmp(varargin{1},'L-SIG')
        ind = indLSIG;
        return;
    end

    nFieldSamples = trc.TRLSIG*sf;
    indRLSIG = uint32([numCumSamples+1 numCumSamples+nFieldSamples]);
    numCumSamples = numCumSamples+nFieldSamples;
    if nargin>1 && strcmp(varargin{1},'RL-SIG')
        ind = indRLSIG;
        return;
    end

    % HE-SIG-A
    nFieldSamples = trc.THESIGA*sf;
    indHESIGA = uint32([numCumSamples+1 numCumSamples+nFieldSamples]);
    numCumSamples = numCumSamples+nFieldSamples;

    if nargin>1 && strcmp(varargin{1},'HE-SIG-A')
        ind = indHESIGA;
        return;
    end

    % HE-SIG-B
    nFieldSamples = 0;
    indHESIGB = zeros(0,2,'uint32');
    numCumSamples = numCumSamples+nFieldSamples;
 
    if nargin>1 && strcmp(varargin{1},'HE-SIG-B')
        ind = indHESIGB;
        return;
    end

    % HE-STF
    nFieldSamples = trc.THESTFNT*sf;
    indHESTF = uint32([numCumSamples+1 numCumSamples+nFieldSamples]);
    numCumSamples = numCumSamples+nFieldSamples;

    if nargin>1 && strcmp(varargin{1},'HE-STF')
        ind = indHESTF;
        return;
    end

    % HE-LTF
    S = heRangingLTFInfo(cfg);
    samplePerSym = trc.THELTFSYM*sf;

    N = [0 cumsum(S.NHELTFWithRepetition)];
    startIndex = numCumSamples+1+N(1:end-1).*samplePerSym;
    endIndex = numCumSamples+N(2:end).*samplePerSym;
    % Append start and end indices of HELTF field
    if cfg.SecureHELTF
        numUsers = numel(cfg.User);
        cpLen = trc.TGIData*sf;
        if numUsers>1
            % The end index includes the zero-power GI after the last user HE-LTF symbol
            indHELTF = uint32([startIndex.' (endIndex+cpLen).']);
        else
            indHELTF = uint32([startIndex endIndex+cpLen]);
        end
    else
        indHELTF = uint32([startIndex endIndex]);
    end
    numCumSamples = endIndex(end);

    if nargin>1 && strcmp(varargin{1},'HE-LTF')
        ind = indHELTF;
        return;
    end

    % HE-Data
    indHEData = zeros(0,2,'uint32');

    % HE-PE
    nFieldSamples = trc.TPE*sf;
    indHEPE = uint32([numCumSamples+1 numCumSamples+nFieldSamples]);

    if nargin>1 && strcmp(varargin{1},'HE-LTF')
        ind = indHELTF;
        return;
    end

    if nargin>1 && strcmp(varargin{1},'HE-Data')
        ind = indHEData;
        return;
    end

    if (nargin>1 && strcmp(varargin{1},'HE-PE'))
        ind = indHEPE;
        return;
    end

    if nargin==1
        % Return indices for all fields
        indField = struct;
        indField.LSTF = indLSTF;
        indField.LLTF = indLLTF;
        indField.LSIG = indLSIG;
        indField.RLSIG = indRLSIG;
        indField.HESIGA = indHESIGA;
        indField.HESIGB = indHESIGB;
        indField.HESTF = indHESTF;
        indField.HELTF = indHELTF;
        indField.HEData = indHEData;
        indField.HEPE = indHEPE;

        ind = indField;
    else
        ind = zeros(0,2,'uint32'); % For codegen
    end
end
