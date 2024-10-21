function positionEst = hePositionEstimate(position,distance)
%hePositionEstimate Estimates the station (STA) position using
%trilateration
%   POSITIONEST = hePositionEstimate(POSITION,DISTANCE) estimates the
%   unknown position of the STA, POSITIONEST using the known access points
%   (APs) positions, POSITION and the estimated distances between APs and
%   STA, DISTANCE.
%
%   POSITION is a matrix of size 2-by-N, represents the position of the N
%   number of APs in a network. Each column of POSITION denotes the
%   2-dimentional position of each AP in xy-plane.
%
%   DISTANCE is a vector of size 1-by-N, represents the distance between the
%   STA and APs in meters.
%
%   POSITIONEST represents the estimated 2-dimentional position of the STA
%   in xy-plane.

%   Copyright 2020 The MathWorks, Inc.

% Load the position of APs in X and Y coordinates
numAPs = size(position,2); % Number of APs
dimLen = size(position,1); % Dimensions (2D/3D)
x = position(1,:).';
y = position(2,:).';

% Perform trilateration to estimate the position of the STA in meters
A = zeros(numAPs-1,dimLen);
B = zeros(1,numAPs-1);
for i = 1:numAPs-1
    A(i,:) = [2*(x(i)-x(end)) 2*(y(i)-y(end))];
    B(i) = x(i)^2-x(end)^2+y(i)^2-y(end)^2+distance(end)^2-distance(i)^2;
end
positionEst = (A.'*A)\(A.'*B(:));
end