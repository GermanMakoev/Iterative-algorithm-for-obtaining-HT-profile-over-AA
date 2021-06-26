filename = 'dip_140_05_16e8.mat';

gridstep = [0.5 0.5]; % arcsec

maxB = 3000; % G

H = [1, 1.5e8, 1.8e8, 2e10]; % cm
Temp = [1e4, 1e4, 2.5e6, 2.5e6]; % K
Dens = [1e12 1e12 1e9 1e9]; % cm^{-3}

freq = 8e9; % Hz

harms = 2:4;

[pFRight, pScanRight, pFLeft, pScanLeft, base, rc] = mfoGCalcMapSimpleCore(filename, H, Temp, Dens, freq, harms, gridstep, maxB);

disp('*** calculated ***');

plot(pScanRight.*10^(4));
hold on
plot(pScanLeft);

