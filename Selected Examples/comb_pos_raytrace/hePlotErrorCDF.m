function hePlotErrorCDF(error,snrRange)
%hePlotErrorCDF Plots the CDF of errors
%
%   hePlotErrorCDF(ERROR,SNRRANGE) plots the empirical CDF of errors for
%   multiple SNR points
%
%   ERROR is a N-by-P array of errors where N is the number of error
%   measurements at a particular SNR and P is the number of SNR points
%
%   SNRRANGE is a P-by-1 array of SNR values for which ERROR was generated

%   Copyright 2020 The MathWorks, Inc.

figure;
totName = [];
xRange = zeros(1,length(snrRange));
for isnr = 1:length(snrRange)
    [x, F] = stairs(sort(error(:,isnr)),(1:length(error(:,isnr)))/length(error(:,isnr)));
    plot(x,F)
    hold on
    [~,idx] = min(abs(F - 0.9));
    xRange(isnr) = x(idx);
    name = ['SNR: ',num2str(snrRange(isnr)),' dB'];
    totName = [convertCharsToStrings(totName),convertCharsToStrings(name)];
end
if isscalar(xRange) && isnan(xRange)
    xlim([0 1]),ylim([0 1]);
else
    xlim([0 max(xRange,[],'omitnan')]),ylim([0 1]);
end
legend(totName,'Location','Southeast');

xlabel('RMS positioning error (meters)'),ylabel('CDF');
title('Positioning Error CDF');
grid on;

end