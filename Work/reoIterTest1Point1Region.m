function reoIterTest1Point1Region

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

gridstep = [1.5 1.5]; % in arcsec
mfoData.posangle = 0;
mfoData.vcos = [0; 0; 1];
M = reoSetField(hLib, mfoData, gridstep);


%--- model ----------------------------------------------------------
freqs = (3:0.3:18)*1e9;
NT = 6e15;

H0 = [1, 1.4e8, 1.6e8, 5e8, 2e10];
Temp0 = [1e4, 1e4, 1e6, 2e6, 2e6];
Dens0 = NT./Temp0;
save('H0.mat','H0');
save('Temp0.mat','Temp0');
Robs = zeros(1,length(freqs));
Lobs = zeros(1,length(freqs));
posR = zeros(1,length(freqs));
posL = zeros(1,length(freqs));
for k = 1:length(freqs)
    [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
     depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft, ...
     pScanRight, pScanLeft] = ...
             reoCalculate(hLib, mfoData, H0, Temp0, Dens0, M, freqs(k), 2:4, 25, 3);
    [Robs(k), posR(k)] = max(pScanRight);
    [Lobs(k), posL(k)] = max(pScanLeft);
    
end

%freqs(51:end) = [];
%Robs(51:end) = [];
%Lobs(51:end) = [];
%--------------------------------------------------------------------------

Ht1 = [1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15]*1e8;
Ht2 = [1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 15 20]*1e8;
Hc = (Ht1+Ht2)/2;
%Tc  = [1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e4 1e6 1e6 1e6 1e6 1e6 1e6 1e6 1e6 1e6];
Tc = 1e6*ones(1, length(Hc));
Tb = 10000;
Hb = 1.2e8;

%Увеличение поля
%FF=0.75;
FF = 0.2;
mfoData.B.x=mfoData.B.x*(1+FF);
mfoData.B.y=mfoData.B.y*(1+FF);
mfoData.B.z=mfoData.B.z*(1+FF);
M = reoSetField(hLib, mfoData, gridstep);

param = reoGetParam;
%param.wTemp = 0.2;
param.wTemp = 20;
param.rescntmax = 30;
%NT = 6e15;
reoIterationCore1Point1Region(hLib, mfoData, M, freqs, Robs, Lobs, H0, Temp0, Hb, Tb, Ht1, Ht2, Hc, Tc, NT, posR, posL, param);

utilsFreeLibrary(hLib);

end
