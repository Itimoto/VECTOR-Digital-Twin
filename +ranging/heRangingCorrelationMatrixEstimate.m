function R = heRangingCorrelationMatrixEstimate(CFR)
%heRangingCorrelationMatrixEstimate Estimates the correlation matrix
%
%   R = heRangingCorrelationMatrixEstimate(CFR) estimates a
%   spatially-smoothed and forward-backward averaged correlation matrix R
%   from a CFR (Channel Frequency Response)
%
%   CFR is a complex matrix of size N-by-P where N is the length of a FFT
%   of a single CFR snapshot and P is the number of snapshots. CFR
%   represents the channel frequency response.
%
%   R is a complex L-by-L matrix representing the estimated correlation
%   matrix.

%   Copyright 2020 The MathWorks, Inc.

% Correlation matrix estimated as per Xinrong Li and K. Pahlavan
% Super-resolution TOA estimation with diversity for indoor geolocation.
% IEEE Transactions on Wireless Communications, vol. 3, no. 1, pp. 224-234,
% Jan. 2004

[N, P] = size(CFR); % P is the number of snapshots, N is the FFT length
L = round(2/3*N); % length of the subarrays for spatial smoothing
M = N-L+1; % number of subarrays
R = complex(zeros(L,L)); % Correlation matrix

if P==1 % Single CFR snapshot
    
    % Perform spatial-smoothing
    for k = 1:M
        subArray = CFR(k:k+L-1);
        R = subArray*subArray' + R;
    end
    R = R/M;
    
else % Multiple CFR snapshots
    
    Rn = complex(zeros(N,N));
    for k = 1:P
        Rn = CFR(:,k)*CFR(:,k)' + Rn;
    end
    Rn = Rn/P;
    
    % Spatial smoothing
    for k = 1:M
        R = R+ Rn(k:k+L-1,k:k+L-1);
    end
    R = R/M;
    
end
J = fliplr(eye(L));
R = (R+J*conj(R)*J)/2; % Forward-backward averaging

end