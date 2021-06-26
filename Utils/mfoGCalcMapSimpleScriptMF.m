filename = 'dip_140_05_16e8.mat';

gridstep = [0.5 0.5]; % in arcsec

maxB = 3230; % G

H = [1, 1.5e8, 1.8e8, 2e10]; % cm
Temp = [1e4, 1e4, 2.5e6, 2.5e6]; % K
Dens = [1e12 1e12 1e9 1e9]; % cm^{-3}

freqs = [3:18]*1e9;

harms = 2:4;

maxR = zeros(1, length(freqs));
maxL = zeros(1, length(freqs));
for k = 1:length(freqs)
    [~, pScanRight, ~, pScanLeft, ~, rc] = mfoGCalcMapSimpleCore(filename, H, Temp, Dens, freqs(k), harms, gridstep, maxB);
    maxR(k) = max(pScanRight);
    maxL(k) = max(pScanLeft);
end

disp('*** calculated ***')
plot(freqs, maxR);
hold on
plot(freqs, maxL);