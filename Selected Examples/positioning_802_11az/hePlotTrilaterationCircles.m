function hePlotTrilaterationCircles(positionAP,positionSTAEst,distanceEst,snr,iter)
%hePlotTrilaterationCircles Plots the positioning trilateration circles
%
%   hePlotTrilaterationCircles(POSITIONAP,POSITIONSTAEST,DISTANCEEST,SNR,ITER)
%   plots the STA, estimated STA and APs positions in the xy-plane with
%   trilateration circles at the given signal-to-noise ratio.
%
%   POSITIONAP represents the 2-dimentional position of the APs in
%   xy-plane.
%
%   POSITIONSTAEST represents the estimated 2-dimentional position of the
%   STA in xy-plane.
%
%   DISTANCEEST represents the estimated distance between the STA and APs
%   in meters.
%
%   SNR is the signal-to-noise ratio in dB.
%
%   ITER denotes the iteration.

%   Copyright 2020 The MathWorks, Inc.

% Plot the positions of the STA and APs in the xy-plane
positionSTA = [0; 0]; %STA always at the origin

figure
hePlotNodePositions(positionAP,positionSTA); % Plot AP and actual STA location
plot(positionSTAEst(1),positionSTAEst(2),'kx','LineWidth',1.5,'MarkerSize',10); % Plot STA estimate

% Plot trilateration circles
numAPs = size(positionAP,2);
angles = 0:2*pi/720:2*pi;
for i = 1:numAPs
    x = distanceEst(i) * cos(angles) + positionAP(1,i);
    y = distanceEst(i) * sin(angles) + positionAP(2,i);
    plot(x,y,'k--','LineWidth',1);
    hold on; grid on
end
title(['Node positions at SNR ' num2str(snr) ' dB for iteration #' num2str(iter) ])
legend({'AP position','STA position','Estimated STA position','Trilateration circles'},'Location','best','FontSize',10)
legend('boxoff')
end