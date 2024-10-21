function [nest,sigest] = ehtNoiseEstimate(x,chanEstSSPilots,cfg,varargin)
%ehtNoiseEstimate Estimate noise power using EHT data field pilots
%
%   NEST = ehtNoiseEstimate(x,CHANESTSSPILOTS,CFG) estimates the mean noise
%   power in watts using the demodulated pilot symbols in the EHT data
%   field and single-stream channel estimates at pilot subcarriers. The
%   noise estimate is averaged over the number of symbols and receive
%   antennas.
%
%   X is a complex Nsp-by-Nsym-by-Nr array containing demodulated pilot
%   subcarrier in EHT data field. Nsym is the number of demodulated
%   EHT-Data symbols.
%
%   CHANESTSSPILOTS is a complex Nsp-by-Nltf-by-Nr array containing the
%   channel gains at pilot subcarrier locations for each symbol, assuming
%   one space-time stream at the transmitter. Nltf is the number of EHT-LTF
%   symbols.
%
%   CFG is a format configuration object of type <a href="matlab:help('wlanEHTMUConfig')">wlanEHTMUConfig</a>,
%   <a href="matlab:help('wlanEHTTBConfig')">wlanEHTTBConfig</a>, or <a href="matlab:help('ehtTBSystemConfig')">ehtTBSystemConfig</a>.
%
%   NEST = ehtNoiseEstimate(...,RUNUMBER) performs noise power estimation
%   for an EHT packet format.
%
%   #  For an EHT MU OFDMA PPDU type, RUNUMBER is required.
%   #  For an EHT MU non-OFDMA PPDU type, RUNUMBER is not required.
%   #  For an EHT TB PPDU type RUNUMBER is not required.
%
%   [NEST,SIGEST] = ehtNoiseEstimate(...) additionally returns an estimate
%   of the signal power.

%   Copyright 2020-2022 The MathWorks, Inc.

narginchk(3,4)
validateattributes(cfg,{'wlanEHTMUConfig','wlanEHTTBConfig','ehtTBSystemConfig'},{'scalar'},mfilename,'format configuration object');

numOFDMSym = size(x,2);
n = 0:numOFDMSym-1;

ruNum = 1;
if isa(cfg,'wlanEHTMUConfig')
    if nargin==4
        ruNum = varargin{1};
    end
    allocInfo = ruInfo(cfg);
    ruSize = allocInfo.RUSizes{ruNum};
    sigInfo = wlan.internal.ehtSIGCodingInfo(cfg);
    numEHTSIG = sigInfo.NumSIGSymbols; % Number of OFDM symbols in EHT-SIG field
else % EHT TB and EHT TB system
    allocInfo = ruInfo(cfg);
    ruSize = allocInfo.RUSizes{ruNum};
    numEHTSIG = 0;
end

numUSIG = 2; % Number of OFDM symbols in U-SIG field

z = 2+numUSIG+numEHTSIG; % Pilot symbol offset
% Get the reference pilots for one space-time stream, pilot sequence same
% for all space-time streams
refPilots = wlan.internal.ehtPilots(ruSize,1,n,z); % Nsp-by-Nsym-by-1

% Average single-stream pilot estimates over symbols (2nd dimension)
avChanEstSSPilots = mean(chanEstSSPilots,2); % Nsp-by-1-by-Nrx

% Estimate channel at pilot location using least square estimates
chanEstPilotsLoc = x./refPilots; % Nsp-by-Nsym-by-Nrx

% Subtract the noisy least squares estimates of the channel at pilot symbol
% locations from the noise averaged single stream pilot symbol estimates of
% the channel
error = chanEstPilotsLoc-avChanEstSSPilots; % Nsp-by-Nsym-by-Nrx

% Get power of error and average over pilot symbols, subcarriers and
% receive antennas
useIdx = ~isnan(error); % NaNs may exist in 1xEHTLTF
nest = real(mean(error(useIdx).*conj(error(useIdx)),'all'));

if nargout>1
   % Get power of channel estimate at pilot locations
   sigest = real(mean(chanEstPilotsLoc(:).*conj(chanEstPilotsLoc(:))));
end
end