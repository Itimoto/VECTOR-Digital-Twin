function [positionSTA,positionAP,radius] = heGeneratePositions(n)
%heGeneratePositions Generate STA and APs position
%
%   [POSITIONSTA,POSITIONAP,RADIUS] = heGeneratePositions(N) places a
%   station (STA) and N number of access points (APs) randomly in xy-plane
%   with X and Y positions in meters.
%
%   POSITIONSTA represents the 2-dimentional position of the STA at origin,
%   with coordinates as [0; 0].
%
%   POSITIONAP is a matrix of size 2-by-N, represents the position of the
%   APs. Each column of POSITIONAP denotes the 2-dimentional position of
%   each AP in xy-plane.
%
%   RADIUS is a vector of size 1-by-N, represents the distance between the
%   STA and APs in meters.
%
%   N is the number of APs in a network.

%   Copyright 2020 The MathWorks, Inc.

% Position the STA
positionSTA = [0; 0]; % STA is always assumed to be at origin

% Generate angle values in radians within [0: 2*pi] with a 2*pi/n radian
% sector
phi = (0:n-1)*(2*pi/n)+rand(1,n)*2*pi/(2*n);

% Generate radius values in meters from the uniform distribution on the
% interval [0,100]
radius = 100.*rand(1,n); % Radius in meters

% Transform polar to Cartesian coordinates
[x, y] = pol2cart(phi,radius);

% Position of APs
positionAP = [x; y]; % X and Y positions in meters

end