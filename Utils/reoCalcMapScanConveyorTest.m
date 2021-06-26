function reoCalcMapScanConveyorTest

mfoData = iouLoadField('C:\Users\Alan Makoev\Desktop\Matlab\dip_140_05_16e8.mat');

if isempty(mfoData)
    return
end

QT = 1;
freefree = 0;

hLib = reoInitLibrary(QT, freefree);
if hLib == 0
    return
end

gridstep = [1 1]; % in arcsec
posangle = 0;
vcos = [0; 0; 1];

[M, base] = reoSetField(hLib, mfoData, gridstep, posangle, vcos);

H = [1, 1e8, 1.1e8, 2e10];
Temp = [1e4, 1e4, 5e6, 5e6];
Dens = 3e15./Temp;
freq = 10.9e9;
harms = 2:4;
pTauL = [1 25];
mode = 2;
c = 0;
b = 0;

[depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
 depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft, ...
 pScanRight, pScanLeft] = ...
         reoCalculate(hLib, mfoData, H, Temp, Dens, M, freq, harms, pTauL, mode, c, b);

[maxR, posR] = max(pScanRight);
[maxL, posL] = max(pScanLeft);
     
surf(pFRight);

end
