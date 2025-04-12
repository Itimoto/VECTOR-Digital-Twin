function [lat2, lon2] = offsetLatLon(lat1, lon1, d, theta)
% Computes the destination point given a start point,
% distance, and bearing using Vincenty's direct formula on the WGS84 ellipsoid.
%
%   [lat2, lon2] = offsetLatLonEllipsoid(lat1, lon1, s, alpha1)
%
%   Inputs:
%       lat1   - initial latitude in degrees
%       lon1   - initial longitude in degrees
%       d      - distance to travel (meters)
%       theta  - initial bearing in degrees (clockwise from north)
%
%   Outputs:
%       lat2   - destination latitude in degrees
%       lon2   - destination longitude in degrees
%
%   The function uses WGS84 ellipsoid parameters:
%       a = 6378137.0 (semi-major axis)
%       f = 1/298.257223563 (flattening)
%
%   Reference:
%       Vincenty, Thaddeus (1975). "Direct and Inverse Solutions of Geodesics on
%       the Ellipsoid with application of nested equations". Survey Review.
%
%   Example:
%       [newLat, newLon] = offsetLatLonEllipsoid(40.7128, -74.0060, 1000, 72);

    % WGS84 ellipsoid parameters
    a = 6378137.0;                 % semi-major axis in meters
    f = 1/298.257223563;           % flattening
    b = (1 - f) * a;               % semi-minor axis
    s = d;

    % Convert input angles from degrees to radians
    phi1 = deg2rad(lat1);
    lambda1 = deg2rad(lon1);
    alpha1_rad = deg2rad(theta);
    
    % U1: 'reduced latitude'
    U1 = atan((1 - f) * tan(phi1));
    
    % Precompute values
    sinU1 = sin(U1);
    cosU1 = cos(U1);
    sinAlpha1 = sin(alpha1_rad);
    cosAlpha1 = cos(alpha1_rad);
    
    % Calculate sigma1, the angular distance on the sphere from the equator to phi1
    sigma1 = atan2( tan(U1), cosAlpha1 );
    
    % Compute sine of alpha and its square
    sinAlpha = cosU1 * sinAlpha1;
    cosSqAlpha = 1 - sinAlpha^2;
    
    % Compute uSq: square of the ellipse parameter
    uSq = cosSqAlpha * (a^2 - b^2) / b^2;
    
    % Coefficients A and B for the series expansion
    A = 1 + (uSq / 16384) * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
    B = (uSq / 1024) * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
    
    % Initial sigma (angular distance on the sphere)
    sigma = s / (b * A);
    
    % Iterate until convergence
    sigmaP = 2 * pi;
    while abs(sigma - sigmaP) > 1e-12
        sigmaM2 = 2 * sigma1 + sigma;  % 2*sigma1 + sigma
        % Compute deltaSigma using the series expansion
        deltaSigma = B * sin(sigma) * ( cos(sigmaM2) + (B/4) * ( cos(sigma) * (-1 + 2*cos(sigmaM2)^2) - ...
                        (B/6) * cos(sigmaM2) * (-3 + 4*sin(sigma)^2) * (-3 + 4*cos(sigmaM2)^2) ) );
        sigmaP = sigma;
        sigma = s / (b * A) + deltaSigma;
    end
    
    % Compute destination latitude
    tmp = sinU1 * cos(sigma) + cosU1 * sin(sigma) * cosAlpha1;
    phi2 = atan2( tmp, (1 - f) * sqrt(sinAlpha^2 + (sinU1 * sin(sigma) - cosU1 * cos(sigma) * cosAlpha1)^2) );
    
    % Compute destination longitude
    lambda = atan2( sin(sigma) * sinAlpha1, cosU1 * cos(sigma) - sinU1 * sin(sigma) * cosAlpha1 );
    % Correction term C
    C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
    L = lambda - (1 - C) * f * sinAlpha * ( sigma + C * sin(sigma) * ( cos(2 * sigma1 + sigma) + C * cos(sigma) * (-1 + 2 * cos(2 * sigma1 + sigma)^2) ) );
    
    % Compute final longitude
    lambda2 = lambda1 + L;
    
    % Convert results back to degrees
    lat2 = rad2deg(phi2);
    lon2 = rad2deg(lambda2);
end