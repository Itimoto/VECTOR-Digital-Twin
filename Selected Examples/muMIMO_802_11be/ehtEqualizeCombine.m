function [y,csi] = ehtEqualizeCombine(x,chanEst,nVar,cfg,varargin)
%ehtEqualizeCombine EHT MIMO channel equalization
%
%   [Y,CSI] = ehtEqualizeCombine(X,CHANEST,NOISEVAR,CFG) performs
%   minimum-mean-square-error (MMSE) frequency domain equalization using
%   the signal input X, the channel estimate, CHANEST, and noise variance,
%   NVAR.
%
%   Y is an estimate of the transmitted frequency domain signal and is of
%   size Nsd-by-Nsym-by-Nsts, where Nsd represents the number of carriers
%   (frequency domain), Nsym represents the number of symbols (time
%   domain), and Nsts represents the number of space-time streams (spatial
%   domain). It is complex when either X or CHANEST is complex, or is real
%   otherwise.
%
%   CSI is a real matrix of size Nsd-by-Nsts containing the soft channel
%   state information.
%
%   X is a real or complex array containing the frequency domain signal to
%   equalize. It is of size Nsd-by-Nsym-by-Nr, where Nr represents the
%   number of receive antennas.
%
%   CHANEST is a real or complex array containing the channel estimates for
%   each carrier and symbol. It is of size Nsd-by-Nsts-by-Nr.
%
%   NVAR is a nonnegative scalar representing the noise variance.
%
%   CFG is a format configuration object of type <a href="matlab:help('wlanEHTMUConfig')">wlanEHTMUConfig</a>,
%   <a href="matlab:help('wlanEHTTBConfig')">wlanEHTTBConfig</a>, or <a href="matlab:help('ehtTBSystemConfig')">ehtTBSystemConfig</a>.
%
%   [Y,CSI] = ehtEqualizeCombine(..., USERIDX) performs noise power
%   estimation for an EHT packet format.
%
%   #  For an EHT MU OFDMA and MU-MIMO, non-OFDMA PPDU type, USERIDX is the
%      1-based index of the user to decode within the EHT MU transmission.
%   #  For an EHT MU, single user, non-OFDMA PPDU type, USERIDX is not
%      required.
%   #  For an EHT TB PPDU type USERIDX is not required.

%   Copyright 2020-2022 The MathWorks, Inc.

validateattributes(cfg,{'wlanEHTMUConfig','wlanEHTTBConfig','ehtTBSystemConfig'},{'scalar'},mfilename,'format configuration object');

% Equalize 
[y,csi] = ofdmEqualize(x,chanEst,nVar);
if isa(cfg,'wlanEHTMUConfig')
    userIdx = 1; % Default
    if nargin>4
        userIdx = varargin{1};
    end
    % Get the indices of the space-time streams for this user
    allSTSIdx = wlan.internal.heSpaceTimeStreamIndices(cfg);
    stsIdx = allSTSIdx(1,userIdx):allSTSIdx(2,userIdx);
    % Extract used STS
    y = y(:,:,stsIdx);
    csi = csi(:,stsIdx);
end