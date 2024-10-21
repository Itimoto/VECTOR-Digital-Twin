function [chanEst,pktOffset] = heRangingSynchronize(rx,cfg)
%heRangingSynchronize Performs Ranging synchronization
%
%   [CHANEST,PKTOFFSET] = heRangingSynchronize(RX,CFG) returns the
%   estimated channel, CHANEST, between all space-time streams and receive
%   antennas using HE-LTF of an HE ranging NDP and PKTOFFSET, an integer
%   scalar indicating the location of the start of a detected packet. If no
%   packet is detected an empty value is returned for both CHANEST and
%   PKTOFFSET.
%
%   RX is the received time-domain signal on which ranging synchronization
%   is performed. It is an Ns-by-Nr matrix of real or complex values, where
%   Ns represents the number of time-domain samples and Nr represents the
%   number of receive antennas.
%
%   CFG is a format configuration object of type <a
%   href="matlab:help('heRangingConfig')">heRangingConfig</a>.
%
%   CHANEST is a complex Nst-by-Nsts-by-Nr array characterizing the
%   estimated channel, where Nst is the number of occupied subcarriers,
%   Nsts is the total number of space-time streams, and Nr is the number of
%   receive antennas.
%
%   PKTOFFSET is an integer scalar indicating the location of the start
%   of a detected packet as the offset from the start of the signal RX. If
%   no packet is detected an empty value is returned.

%   Copyright 2020-2021 The MathWorks, Inc.

% Get channel bandwidth
chanBW = cfg.ChannelBandwidth;

% Get sample rate and the ranging field indices
fs = wlanSampleRate(chanBW);
ind = ranging.heRangingFieldIndices(cfg);

% Packet detect and determine coarse packet offset
coarsePktOffset = wlanPacketDetect(rx,chanBW);
if ~isempty(coarsePktOffset) && ((coarsePktOffset+ind.LSIG(2))<=size(rx,1))
    % Extract L-STF and perform coarse frequency offset correction
    lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
    coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
    rx = frequencyOffset(rx,fs,-coarseFreqOff);  
    
    % Extract the non-HT fields and determine fine packet offset
    nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
    finePktOffset = wlanSymbolTimingEstimate(nonhtfields,chanBW);

    % Determine final packet offset
    pktOffset = coarsePktOffset+finePktOffset;
    
    % Check if the packet is detected out of the range. With a maximum
    % propagation distance of 100 meters, the expected maximum propagation
    % delay should be 54 samples for channel bandwidth of 160 MHz.
    % maxPktOffset = (160e6*100 meters/speedOfLight)
    if (pktOffset <= 54) && ((pktOffset+ind.HELTF(2))<=size(rx,1))
        % Extract L-LTF and perform fine frequency offset correction
        rxLLTF = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
        fineFreqOff = wlanFineCFOEstimate(rxLLTF,chanBW);
        rx = frequencyOffset(rx,fs,-fineFreqOff);
        rxHELTF = rx(pktOffset+(ind.HELTF(1):ind.HELTF(2)),:);

        % OFDM demodulate
        heltf = ranging.heRangingDemodulate(rxHELTF,cfg);

        % Channel estimate
        chanEst = ranging.heLTFRangingChannelEstimate(heltf,cfg);
        
    else
        chanEst = [];
        pktOffset = [];
    end
else
    chanEst = [];
    pktOffset = [];    
end
end
