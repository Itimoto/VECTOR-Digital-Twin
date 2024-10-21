function heRangingPlotPositions(positionAP,positionSTAEst,distanceEst,snr,iter)
%heRangingPlotPositions Plots the HE ranging positioning results
%   heRangingPlotPositions(POSITIONAP,POSITIONSTAEST,DISTANCEEST,SNR,ITER)
%   plots the STA, estimated STA and APs positions in the xy-plane with
%   trilateration circles at the given signal-to-noise ratio.
%
%   POSITIONAP represents the 2-dimentional position of the APs in xy-plane.
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

positionSTA = [0;0];
numAPs = size(positionAP,2);
angles = 0:2*pi/720:2*pi;
figure
for i = 1:numAPs
    x = distanceEst(i) * cos(angles) + positionAP(1,i);
    y = distanceEst(i) * sin(angles) + positionAP(2,i);
    h(1) = plot(x,y,'k:','LineWidth',2);
    hold on; grid on
    h(2) = plot(positionAP(1,i),positionAP(2,i),'b*',...
        'MarkerSize',5,'LineWidth',1);
end
h(3) = plot(positionSTA(1),positionSTA(2),'rP','MarkerSize',10,'LineWidth',2);
h(4) = plot(positionSTAEst(1),positionSTAEst(2),'kS','MarkerSize',10,'LineWidth',2);
xlabel('X-position (meters)')
ylabel('Y-position (meters)')
title(['Node positions at SNR ' num2str(snr) ' dB for iteration #' num2str(iter) ])
legend(h([1 2,3,4]),{'Trilateration circles','AP position','STA position','Estimated STA position'},'Location','best','FontSize',10)
legend('boxoff')
axis equal
end