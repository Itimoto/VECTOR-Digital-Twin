classdef heRangingUser < comm.internal.ConfigBase
%heRangingUser User properties of HE-LTF field
%   CFGUser = heRangingUser creates a user configuration object. This
%   object contains the HE-LTF properties of a user within an HE Ranging
%   NDP PPDU.
%
%   CFGUSER = heRangingUser(Name,Value) creates an object that holds
%   the HE-LTF properties for the users within an HE Ranging NDP PPDU,
%   CFGUSER, with the specified property Name set to the specified value.
%   You can specify additional name-value pair arguments in any order as
%   (Name1,Value1, ...,NameN,ValueN).
%
%   heRangingUser objects are used to parameterize users within an HE
%   Ranging NDP PPDU transmission, and therefore are part of the
%   <a href="matlab:help('heRangingConfig')">heRangingConfig</a> object.
%
%   heRangingUser properties:
%
%   NumSpaceTimeStreams - Number of space-time streams
%   NumHELTFRepetition  - Number of repetitions of HE-LTF symbols
%   SecureHELTFSequence - Secure HE-LTF sequence
%
%   See also heRangingConfig.

%   Copyright 2020-2022 The MathWorks, Inc.

%#codegen

properties
    %NumSpaceTimeStreams Number of space-time streams
    %   Specify the number of space-time streams as integer between 1 and
    %   8, inclusive. The default value of this property is 1.
    NumSpaceTimeStreams (1,1) {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(NumSpaceTimeStreams,1),mustBeLessThanOrEqual(NumSpaceTimeStreams,8)} = 1;
    %NumHELTFRepetition Number of repetitions of HE-LTF symbols
    %   Specify the number of HE-LTF repetition as integer between 1 and 8,
    %   inclusive. The default value of this property is 1.
    NumHELTFRepetition (1,1) {mustBeNumeric,mustBeInteger,mustBeGreaterThanOrEqual(NumHELTFRepetition,1),mustBeLessThanOrEqual(NumHELTFRepetition,8)} = 1;
    %SecureHELTFSequence Secure HE-LTF sequence
    %   Specify the secure HE-LTF sequence as a character vector or string
    %   scalar representing octets in hexadecimal format. The maximum
    %   number of bits corresponding to the hexadecimal sequence are
    %   derived from Equation 11-yy and 11-zz of IEEE P802.11az/D2.0,
    %   Section 11.22.6.4.6.3. The sequence is specified for a single
    %   measurement instance. If the number of input bits in the
    %   hexadecimal sequence is less than the required number of bits for
    %   the given user configuration, the secure sequence is cyclically
    %   extended. If the number of input bits in the hexadecimal sequence
    %   is more than the required number of bits for the given user
    %   configuration, only the required number of bits are extracted from
    %   the input hexadecimal sequence. The <a href="matlab:help('heRangingConfig/numSecureHELTFBits')">numSecureHELTFBits</a> method of
    %   <a href="matlab:help('heRangingConfig')">heRangingConfig</a> object returns the number of secure HE-LTF bits 
    %   required to generate the secure HE-LTF symbols for a given
    %   configuration. The default sequence is '00000000'.
    SecureHELTFSequence = '00000000';
end

methods
    function obj = heRangingUser(varargin)
        % For codegen set different dimensions to make sure varsize
        if ~isempty(coder.target)
            secureHELTFSequence = '00000000';
            coder.varsize('SecureHELTFSequence',[1 344],[0 1]); % Add variable-size support
            obj.SecureHELTFSequence = secureHELTFSequence; % Default
        end
        obj = setProperties(obj,varargin{:}); % Supperclass method for NV pair parsing
    end

    function obj = set.SecureHELTFSequence(obj,val)
        propName = 'SecureHELTFSequence';

        validateattributes(val,{'char','string'},{},mfilename,propName);
        val = char(upper(val));
        % Validate hex digits
        coder.internal.errorIf(any(~(((val>='0') & (val<='9')) | ((val>='A') & (val<= 'F')))), ...
            'wlan:shared:InvalidHexDigit',propName);

        % Length of hexadecimal data octets must be multiple of 2
        coder.internal.errorIf((rem(numel(val),2)~=0),'wlan:shared:HexNibbleMissing',val);
        obj.(propName) = val;
    end
end
end
