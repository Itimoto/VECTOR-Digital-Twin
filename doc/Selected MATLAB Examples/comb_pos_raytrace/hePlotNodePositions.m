function hePlotNodePositions(positionAP,positionSTA)
%hePlotNodePositions Plots the position of all APs and the STA
%
%   hePlotNodePositions(POSITIONAP,POSITIONSTA) plots the STA and all AP
%   positions in the xy-plane.
%
%   POSITIONAP represents the 2-dimentional position of the APs in
%   xy-plane.
%
%   POSITIONSTA represents the 2-dimentional position of the STA in
%   xy-plane.

%   Copyright 2020 The MathWorks, Inc.

plot(positionAP(1,:),positionAP(2,:),'b^','LineWidth',2),hold on
plot(positionSTA(1),positionSTA(2),'rp','LineWidth',2),hold on
legend('AP Locations','STA Location')

% Get distance from APs to STA to keep consistent plot axes
radiusAP = sqrt(positionAP(1,:).^2 + positionAP(2,:).^2);
range = max(abs(positionAP+repmat(radiusAP,[2 1])),[],'all');
xlim([-range range]),ylim([-range range]);
axis equal
grid on;
xlabel('x-axis (meters)'),ylabel('y-axis (meters)');

end