classdef heRangingConfig < comm.internal.ConfigBase
%heRangingConfig Creates an object for parameterizing HE ranging NDP transmission
%   CFG = heRangingConfig creates an object for parameterizing HE ranging
%   null data packet (NDP) transmission for an HE NDP packet. This object
%   contains transmit parameters for the HE Ranging NDP format of IEEE
%   P802.11az/D2.0 standard.
%
%   CFG = heRangingConfig(NumUsers) creates an object for parameterizing HE
%   ranging NDP transmission for an HE NDP packet. NumUsers specifies the
%   number of users as scalar between 1 and 64 (inclusive). When NumUsers
%   is greater than one the HE-LTF field contains secure HE-LTF symbols for
%   multiple users.
%
%   CFG = heRangingConfig(..., Name,Value) creates an object for
%   parameterizing HE ranging NDP transmission for an HE NDP packet, CFG,
%   with the specified property Name set to the specified Value. You can
%   specify additional name-value pair arguments in any order as
%   (Name1,Value1, ...,NameN,ValueN).
%
%   heRangingConfig methods:
%
%   packetFormat         - HE SU packet format
%   numHELTFSymbols      - Number of HE-LTF symbols
%   numSecureHELTFBits   - Number of secure HE-LTF bits
%   ruInfo               - Resource unit (RU) allocation information
%   showAllocation       - Shows the RU allocation
%   transmitTime         - Time required to transmit a packet
%
%   heRangingConfig properties:
%
%   ChannelBandwidth     - Channel bandwidth
%   InactiveSubchannels  - Indicate punctured 20 MHz subchannels
%   User                 - User properties of HE-LTF
%   SecureHELTF          - Type of HE-LTF sequence
%   NumTransmitAntennas  - Number of transmit antennas
%   PreHECyclicShifts    - Cyclic shift values for >8 transmit chains
%   PreHESpatialMapping  - Indicate spatial mapping of pre-HE-STF portion
%   GuardInterval        - Guard interval type
%   UplinkIndication     - Indicate if the PPDU is sent on the uplink
%   BSSColor             - Basic service set (BSS) color identifier
%   TXOPDuration         - Duration information for TXOP protection
%   HighDoppler          - Indicate high-Doppler mode
%   MidamblePeriodicity  - Midamble periodicity in number of OFDM symbols
%   STBC                 - Enable space-time block coding
%   MCS                  - Modulation and coding scheme
%   DCM                  - Enable dual coded modulation for HE-Data field
%   ChannelCoding        - Channel coding type
%   HELTFType            - Indicate HE-LTF compression mode
%   Beamforming          - Indicate beamforming
%   APEPLength           - Indicate APEP length
%   SpatialReuse         - Indicate spatial reuse indication
%
%   % Example 1:
%   %  Parameterize an HE ranging NDP object for a single user.
%
%   cfg = heRangingConfig;
%   disp(cfg)
%
%   % Example 2:
%   %  Parameterize a secure HE Ranging NDP object for a single user.
%
%   cfg = heRangingConfig;
%   cfg.SecureHELTF = true;
%   disp(cfg)
%
%   % Example 3:
%   %  Parameterize a secure HE Ranging NDP for two users. The number of
%   %  space-time streams are two and one for user 1 and 2, respectively.
%
%   cfg = heRangingConfig(2);    % Create format configuration
%   cfg.NumTransmitAntennas = 2; % Number of transmit antennas
%
%   % Set space-time streams and number of HE-LTF repetitions
%   cfg.User{1}.NumSpaceTimeStreams = 2;
%   cfg.User{1}.NumHELTFRepetition = 2;
%
%   cfg.User{2}.NumSpaceTimeStreams = 1;
%   cfg.User{2}.NumHELTFRepetition = 2;
%
%   % Get the length of the secure sequence in nibbles for both users
%   numNibbles = numSecureHELTFBits(cfg)/4;
%
%   % Set secure HE-LTF sequence for user-1
%   secureSeqUser1 = 'a135c6de1f8a7c3b23b51d2a786a93ab';
%   cfg.User{1}.SecureHELTFSequence = secureSeqUser1(1:numNibbles(1));
%
%   % Set secure HE-LTF sequence for user-2
%   secureSeqUser2 = '17a4b5c256de34fab';
%   cfg.User{2}.SecureHELTFSequence = secureSeqUser2(1:numNibbles(2));
%
%   disp(cfg);
%
%   See also heRangingWaveformGenerator, heRangingFieldIndices

%   Copyright 2020-2022 The MathWorks, Inc.

%#codegen

properties (Access = 'public')
    %ChannelBandwidth Channel bandwidth (MHz) of PPDU transmission
    %   Specify the channel bandwidth as one of 'CBW20' | 'CBW40' | 'CBW80'
    %   | 'CBW160'. The default value of this property is 'CBW20'.
    ChannelBandwidth = 'CBW20';
    %InactiveSubchannels Indicates punctured 20 MHz subchannels
    %   Specify inactive 20 MHz subchannels as a logical vector. The number
    %   of elements must be 1 or equal to the number of 20 MHz subchannels
    %   given ChannelBandwidth. Set an element to true if a 20 MHz
    %   subchannel is inactive (punctured). Subchannels are ordered from
    %   lowest to highest absolute frequency. If a scalar is provided, this
    %   value is assumed for all subchannels. This property applies only
    %   when APEPLength is 0 and ChannelBandwidth is 'CBW80' or 'CBW160'.
    %   At least one subchannel must be active. The default value for this
    %   property is false.
    InactiveSubchannels logical = false;
    %User User properties of HE-LTF
    %   Set the transmission properties of HE-LTF for each user in
    %   the transmission. This property is a cell array of <a
    %   href="matlab:help('heRangingUser')">heRangingUser</a> 
    %   objects. Each element of the cell array contains properties to
    %   configure the HE-LTF for a user.
    User;
    %SecureHELTF Type of HE-LTF sequence
    %   Set this property to true to enable the generation of HE-LTF
    %   symbols with secure HE-LTF sequence as defined in IEEE
    %   P802.11az/D2.0, Section 27.3.17d. When this property is set to
    %   false, regular HE-LTFs as defined in IEEE P802.11ax/D4.1, Section
    %   27.3.10.10 are used. This property is only applicable when the
    %   HE-LTF contains symbols for single user. When the transmission
    %   contains more than one user, the object sets this property to true.
    %   The default is false.
    SecureHELTF (1,1) logical = false;
    %NumTransmitAntennas Number of transmit antennas
    %   Specify the number of transmit antennas as a positive integer
    %   scalar. The default is 1.
    NumTransmitAntennas (1,1) {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(NumTransmitAntennas,1)} = 1;
    %PreHECyclicShifts Cyclic shift values for >8 transmit chains
    %   Specify the cyclic shift values for the pre-HE portion of the
    %   waveform, in nanoseconds, for >8 transmit antennas as a row vector
    %   of length L = NumTransmitAntennas-8. The cyclic shift values must
    %   be between -200 and 0 inclusive. The first 8 antennas use the
    %   cyclic shift values defined in Table 21-10 of IEEE Std 802.11-2016.
    %   The remaining antennas use the cyclic shift values defined in this
    %   property. If the length of this row vector is specified as a value
    %   greater than L the object only uses the first L, PreHECyclicShifts
    %   values. For example, if you specify the NumTransmitAntennas
    %   property as 16 and this property as a row vector of length N>L, the
    %   object only uses the first L = 16-8 = 8 entries. This property
    %   applies only when you set the NumTransmitAntennas property to a
    %   value greater than 8 and PreHESpatialMapping to false. The default
    %   value of this property is -75.
    PreHECyclicShifts {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(PreHECyclicShifts,-200),mustBeLessThanOrEqual(PreHECyclicShifts,0)} = -75;
    %PreHESpatialMapping Indicate spatial mapping of pre-HE-STF portion
    %   Set this property to true to spatially map the pre-HE-STF portion
    %   of the PPDU the same way as the first symbol of the HE-LTF field on
    %   each tone. Set to false to apply no spatial mapping to the
    %   pre-HE-STF portion of the PPDU. The default value of this property
    %   is false.
    PreHESpatialMapping (1,1) logical = false;
    %GuardInterval Guard interval type
    %   Specify the guard interval (cyclic prefix) type in microseconds as
    %   one of 0.8 or 1.6. The default is 1.6.
    GuardInterval {mustBeNumeric,mustBeMember(GuardInterval,[0.8 1.6])} = 1.6;
    %UplinkIndication Indicate uplink transmission
    %   Set this property to true to indicate that the PPDU is sent on an
    %   uplink transmission. The default value of this property is false,
    %   which indicates a downlink transmission.
    UplinkIndication (1,1) logical = false;
    %BSSColor Basic service set (BSS) color identifier
    %   Specify the BSS color identifier of a basic service set as an
    %   integer scalar between 0 to 63, inclusive. The default is 0.
    BSSColor (1,1) {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(BSSColor,0),mustBeLessThanOrEqual(BSSColor,63)} = 0;
    %TXOPDuration Duration information for TXOP protection
    %   Specify the TXOPDuration signaled in HE-SIG-A as an integer scalar
    %   between 0 and 127, inclusive. The TXOP field in HE-SIG-A is set
    %   directly to TXOPDuration, therefore a duration in microseconds must
    %   be converted before being used as specified in Table 27-18 of IEEE
    %   P802.11ax/D4.1. For more information see the <a href="matlab:doc('wlanHESUConfig')">documentation</a>.
    TXOPDuration (1,1) {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(TXOPDuration,0),mustBeLessThanOrEqual(TXOPDuration,127)} = 127;
    %HighDoppler Indicate high-Doppler mode
    %   Set this property to true to indicate high doppler in HE-SIG-A.
    %   This property is only signaled in HE-SIG-A field as HE Ranging
    %   NDP transmission contains no data field. The default value of this
    %   property is false.
    HighDoppler (1,1) logical = false;
    %MidamblePeriodicity Midamble periodicity in number of OFDM symbols
    %   Specify HE-Data field midamble periodicity as 10 or 20 OFDM
    %   symbols. This property applies only when HighDoppler property is
    %   set to true. This property is only signaled in HE-SIG-A field as HE
    %   Ranging NDP transmission contains no data field. The default is 10.
    MidamblePeriodicity (1,1) {mustBeNumeric,mustBeMember(MidamblePeriodicity,[10 20])} = 10;
    %STBC Enable space-time block coding
    %   Set this property to true to enable space-time block coding (STBC)
    %   in the data field transmission. STBC can only be applied for two
    %   space-time streams and when DCM is not used. This property is only
    %   signaled in HE-SIG-A field as HE Ranging NDP transmission contains
    %   no data field. The default value of this property is false.
    STBC (1,1) logical = false;
    %MCS Modulation and coding scheme
    %   Specify the modulation and coding scheme as an integer scalar. Its
    %   elements must be integers between 0 and 11, inclusive. This
    %   property is only signaled in HE-SIG-A field as HE Ranging NDP
    %   transmission contains no data field. The default value of this
    %   property is 0.
    MCS (1,1) {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(MCS,0),mustBeLessThanOrEqual(MCS,11)} = 0;
    %DCM Enable dual-carrier modulation for HE-Data field
    %   Set this property to true to indicate that dual-carrier modulation
    %   (DCM) is used for the HE-Data field. DCM can only be used with MCS
    %   0, 1, 3, 4, when STBC is not used and when the number of space-time
    %   streams is less than or equal to 2. This property is only signaled
    %   in HE-SIG-A field as HE Ranging NDP transmission contains no data
    %   field. The default value of this property is false.
    DCM (1,1) logical = false;
    %ChannelCoding Forward error correction coding type
    %   Specify the channel coding as one of 'BCC' or 'LDPC' to indicate
    %   binary convolution coding (BCC) or low-density-parity-check (LDPC)
    %   coding, respectively. This property is only signaled in HE-SIG-A
    %   field as HE Ranging NDP transmission contains no data field. The
    %   default is 'LDPC'.
    ChannelCoding = 'LDPC';
end

properties (SetAccess = private)
    %HELTFType HE-LTF compression mode of HE Ranging NDP
    %   This is a read-only property. Only 2xHE-LTF compression mode is
    %   supported.
    HELTFType = 2;
    %Beamforming Beamforming indication
    %   This is a read-only property. Beamforming is not supported.
    Beamforming = false;
    %APEPLength APEP length
    %   This is a read-only property. Only APEPLength of 0 is supported.
    APEPLength = 0;
    %SpatialReuse Spatial reuse indication
    %   This is a read-only property. Only SpatialReuse of 15 is
    %   supported.
    SpatialReuse = 15;
end

properties(Dependent,Hidden)
    %SpatialMappingMatrix Spatial mapping matrix
    %   Spatial mapping matrix is of size Nsts-by-Nt, where Nsts
    %   is the maximum number of space-time streams between all users and
    %   Nt is the number of transmit antennas. The SpatialMappingMatrix is
    %   an identity matrix when Nsts=Nt. When Nsts<Nt, the
    %   SpatialMappingMatrix is an identity matrix with all zero rows
    %   removed.
    SpatialMappingMatrix;
    %NumSpaceTimeStreams Maximum number of space-time streams
    %   The maximum number of space-time streams for all users.
    NumSpaceTimeStreams (1,1);
end

properties(Constant,Hidden)
    NominalPacketPadding = 0;
    SpatialMapping = 'Custom';
    ChannelBandwidth_Values = {'CBW20','CBW40','CBW80','CBW160'};
    ChannelCoding_Values = {'BCC','LDPC'};
end

methods
    function obj = heRangingConfig(varargin)
        narginchk(0,Inf);
        if mod(nargin,2)==1 % Odd number of inputs
            numUsers = varargin{1};
            validateattributes(numUsers(1),{'numeric'},{'scalar','>=',1,'<=',64});
            startInd = 2;
        else % Even number of inputs
             numUsers = 1;
             coder.internal.errorIf((mod(nargin,2)~=0),'wlan:ConfigBase:InvalidPVPairs');
             if nargin>0 && ~(ischar(varargin{1})||isstring(varargin{1}))
                error('Invalid input combination');
             end
             startInd = 1;
        end
        for i=startInd:2:nargin-1
             if strcmp(varargin{i},'User')
                error('You cannot set the User property with a Name-Value pair');
             end
        end

        % For codegen set maximum dimensions to force varsize
        if ~isempty(coder.target)
            channelBandwidth = 'CBW20';
            coder.varsize('channelBandwidth',[1 6],[0 1]); % Add variable-size support
            obj.ChannelBandwidth = channelBandwidth; % Default

            channelCoding = 'LDPC';
            coder.varsize('channelCoding',[1 4],[0 1]); % Add variable-size support
            obj.ChannelCoding = channelCoding; % Default
        end
        obj = setProperties(obj,varargin{startInd:end}); % Supperclass method for NV pair parsing

        obj.User = createUsers(numUsers(1)); % For codegen

        % Set SecureHELTF to true when number of users > 1
        if numUsers>1
            obj.SecureHELTF = true;
        end
    end

    function obj = set.ChannelBandwidth(obj,val)
        val = validateEnumProperties(obj,'ChannelBandwidth',val);
        obj.ChannelBandwidth = val;
    end

    function obj = set.ChannelCoding(obj,val)
        val = validateEnumProperties(obj,'ChannelCoding',val);
        obj.ChannelCoding = val;
    end

	function val = get.SpatialMappingMatrix(obj)
        S = ruInfo(obj);
        % Set spatial mapping matrix to an identity matrix as defined in
        % IEEE P802.11az/D2.0, Section 27.3.17a.
        val = eye(S.NumSpaceTimeStreamsPerRU,obj.NumTransmitAntennas);
    end

    function val = get.NumSpaceTimeStreams(obj)
        % Maximum number of space-time stream for all users
        val = maxNumUserSpaceTimeStreams(obj);
    end

    function obj = set.User(obj,value)
    %Validate User property is a cell array of heRangingUser objects

        validateattributes(value,{'cell'},{'nonempty'},mfilename,'User');

        numInp = numel(value);
        for i = 1:numInp
            coder.internal.errorIf(~isa(value{i},'ranging.heRangingUser'),'wlan:shared:InvalidUser',numInp,'ranging.heRangingUser');
        end
        obj.User = value;
    end

    function format = packetFormat(obj) %#ok<MANU>
    %packetFormat Returns the packet format
    %   Returns the packet format as a character vector. Packet format is
    %   'HE-SU'.

       format = 'HE-SU';
    end

    function t = transmitTime(obj,varargin)
    %transmitTime Returns the time required to transmit a packet
    %   T = transmitTime(CFG) returns the time required to transmit a
    %   packet in seconds.
    %
    %   T = transmitTime(CFG,UNIT) returns the transmit time in the
    %   requested unit. UNIT must be 'seconds', 'milliseconds',
    %   'microseconds', or 'nanoseconds'.

        narginchk(1,2);
        [~,t] = wlan.internal.hePLMETxTimePrimative(obj);
        t = t*1e-3; % Scale from nanoseconds to microseconds
        t = wlan.internal.convertTransmitTime(t,varargin{:});
    end

    function numHELTFSym = numHELTFSymbols(obj)
    %numHELTFSymbols Number of HE-LTF symbols
    %   Returns the number of HE-LTF symbols for the given object
    %   configuration

        hezInfo = ranging.heRangingLTFInfo(obj);
        numHELTFSym = hezInfo.NHELTF;
    end

    function numSecureBits = numSecureHELTFBits(obj)
    %numSecureHELTFBits Number of secure HE-LTF bits
    %   Returns the number of secure HE-LTF bits required to generate the
    %   secure HE-LTF symbols for the given configuration. The number of
    %   secure bits are rounded to the nearest multiple of 8.

        hezInfo = ranging.heRangingLTFInfo(obj);
        octetLength = 8; % 1 octet = 8 bits
        numSecureBits = round(hezInfo.NumSecureHELTFBitsPerUser/octetLength)*octetLength;
    end

    function s = ruInfo(obj)
    %ruInfo Returns information relevant to the resource unit
    %   S = ruInfo(cfgHEz) returns a structure, S, containing the resource
    %   unit (RU) allocation information for the heRangingConfig object,
    %   cfgHEz. The output structure S has the following fields:
    %
    %   NumUsers                 - Number of users (1)
    %   NumRUs                   - Number of RUs (1)
    %   RUIndices                - Index of the RU (1)
    %   RUSizes                  - Size of the RU
    %   NumUsersPerRU            - Number of users per RU (1)
    %   NumSpaceTimeStreamsPerRU - Total number of space-time streams
    %   PowerBoostFactorPerRU    - Power boost factor (1)
    %   RUNumbers                - RU number (1)

        s = struct;
        s.NumUsers = 1;
        s.NumRUs = 1;
        s.RUIndices = 1;
        s.RUSizes = wlan.internal.heFullBandRUSize(obj.ChannelBandwidth);
        s.NumUsersPerRU = 1;
        s.NumSpaceTimeStreamsPerRU = maxNumUserSpaceTimeStreams(obj);
        s.PowerBoostFactorPerRU = 1;
        s.RUNumbers = 1;
    end

    function showAllocation(obj,varargin)
    %showAllocation Shows the RU allocation
    %   showAllocation(obj) shows the RU allocation for an HE Ranging NDP
    %   packet format.
    %
    %   showAllocation(obj,AX) shows the allocation in the axes specified
    %   by AX instead of in the current axes.

        wlan.internal.validateInactiveSubchannels(obj);
        wlan.internal.hePlotAllocation(obj,varargin{:});
    end

    function varargout = validateConfig(obj)
    %validateConfig Validate the dependent properties of heRangingConfig object
    %   validateConfig(obj) validates the dependent properties for the
    %   specified heRangingConfig configuration object.
    %
    %   For INTERNAL use only, subject to future changes
    %

        nargoutchk(0,1);

        % Validate PreHECyclicShifts against NumTransmitAntennas
        if ~isInactiveProperty(obj,'PreHECyclicShifts')
            coder.internal.errorIf(numel(obj.PreHECyclicShifts)<obj.NumTransmitAntennas-8, ...
                'wlan:shared:InvalidCyclicShift','PreHECyclicShifts',obj.NumTransmitAntennas-8);
        end

        % Validate inactive sub channels
        if ~isInactiveProperty(obj,'InactiveSubchannels')
            % Puncturing is only applicable for HE SU NDP
            wlan.internal.validateInactiveSubchannels(obj);
        end

        % Validate number of HE-LTF symbols
        if numHELTFSymbols(obj)>64
            error('Number of HE-LTF symbols should be less than or equal to 64');
        end

        % Validate secure HE-LTF symbols against NumUsers. For secure HE-LTF
        % the NumUsers must be greater than 1.
        if numel(obj.User)>1 && ~obj.SecureHELTF
            error('SecureHELTF must be true when NumUser is greater than 1');
        end

        % Validate NumTransmitAntennas and NumSpaceTimeStreams
        if maxNumUserSpaceTimeStreams(obj)>obj.NumTransmitAntennas
        	error('NumSpaceTimeStreams for any user must not be larger than NumTransmitAntennas');
        end

        % Validate MCS and length
        if nargout == 1
            [psduLength,txTime,commonCodingParams] = wlan.internal.hePLMETxTimePrimative(obj);
            sf = 1e3; % Scaling factor to convert time from ns to us
            % Set output structure
            s = struct( ...
                'NumDataSymbols', commonCodingParams.NSYM, ...
                'TxTime', txTime/sf, ...% TxTime in us
                'PSDULength', psduLength);
                varargout{1} = s;
        end
    end
end

methods (Access = protected)
    function flag = isInactiveProperty(obj, prop)
        switch prop
            case 'MidamblePeriodicity'
                % Hide MidamblePeriodicity when HighDoppler is not set
                flag = obj.HighDoppler == false;
            case 'InactiveSubchannels'
                % Hide InactiveSubchannels for non-NDP packet when ChannelBandwidth is CBW20 and CBW40
                flag = ~any(strcmp(obj.ChannelBandwidth,{'CBW80','CBW160'}));
            case 'PreHECyclicShifts'
                % Hide PreHECyclicShifts when NumTransmitAntennas <=8 or pre-HE
                % spatial mapping is used
                flag = obj.NumTransmitAntennas <= 8 || obj.PreHESpatialMapping == true;
            otherwise
                flag = false;
        end
    end
end

end

function user = createUsers(numUsers)
%createUsers User object for HE-LTF field

    % Returns a cell array of Users
    Usertmp = cell(1,numUsers);
    for userIdx = 1:numUsers
        Usertmp{userIdx} = ranging.heRangingUser();
    end
    user = Usertmp;
end

function maxNumSTS = maxNumUserSpaceTimeStreams(cfg)
%maxNumUserSpaceTimeStreams Maximum number of space-time streams for all users

    numUser = numel(cfg.User);
    numSTSPerUser = coder.nullcopy(zeros(1,numUser));
    for u=1:numUser
        numSTSPerUser(u) = cfg.User{u}.NumSpaceTimeStreams;
    end
    % Get the maximum number of space-time streams for all users
    maxNumSTS = max(numSTSPerUser);
end
