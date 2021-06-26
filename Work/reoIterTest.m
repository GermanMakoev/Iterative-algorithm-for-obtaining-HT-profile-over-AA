function reoIterTest

mfoData = iouLoadField('Data\dip_140_05_16e8.mat');
%mfoData = iouLoadField('s:\UData\Data\Dipole\dip_140_05_16e8.mat');
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
% mfoData.vcos = [0; cosd(60); cosd(30)];
mfoData.vcos = [0; 0; 1];
[M, base, Bc] = reoSetField(hLib, mfoData, gridstep);
Mask = reoGetMaskA(sqrt(Bc.x.^2 + Bc.y.^2 + Bc.z.^2));
nModels=2;

%--- model ----------------------------------------------------------
freqs = (3:0.3:18)*1e9;
[diagrH, diagrV] =  mfoCreateRATANDiagrams(freqs, M, gridstep, base, 3);
NT = 3e15;

H0 = [1, 1.2e8, 1.5e8, 1.8e8, 5e8, 2e10];
T0 = [[1e4, 1e4, 1.5e6, 1.8e6, 5e6, 5e6];[1e4, 1e4, 1e4, 2e6, 5e6, 5e6]];
Dens= NT./T0;

Robs = zeros(1,length(freqs));
Lobs = zeros(1,length(freqs));
posR = zeros(1,length(freqs));
posL = zeros(1,length(freqs));
depthRightW=zeros(M);
pFRightW = zeros(M);
pFLRightW = zeros(M);
pHLRightW=zeros(M);
depthLeftW=zeros(M);
pFLeftW = zeros(M);
pFLLeftW = zeros(M);
pHLLeftW=zeros(M);
for k = 1:length(freqs)
            for m = 1:nModels
                 [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
                 depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft] = ...
                         reoCalculate(hLib, mfoData, H0, T0(m,:), Dens(m,:), M, freqs(k), 2:4, 25, 3);
                depthRightW(Mask==m) = depthRight(Mask==m); 
                pFRightW(Mask == m) = pFRight(Mask == m);
                pR = pFLRight(:, :, 1);
                pFLRightW(Mask == m) = pR(Mask == m);
                pHR=pHLRight(:, :, 1);
                pHLRightW(Mask == m) = pHR(Mask == m);
                depthLeftW(Mask==m) = depthLeft(Mask==m); 
                pFLeftW(Mask == m) = pFLeft(Mask == m);
                pL = pFLLeft(:, :, 1);
                pFLLeftW(Mask == m) = pL(Mask == m);
                pHL=pHLLeft(:, :, 1);
                pHLLeftW(Mask == m) = pHL(Mask == m);
            end
            scanR = gstMapConv(pFRightW, diagrH(k,:), diagrV(k,:), gridstep);
            scanL = gstMapConv(pFLeftW, diagrH(k,:), diagrV(k,:), gridstep);
    [Robs(k), posR(k)] = max(scanR);
    [Lobs(k), posL(k)] = max(scanL);
end

%--------------------------------------------------------------------------

Ht1 = [1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15]*1e8;
Ht2 = [1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 15 20]*1e8;
%Ht1 = [1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15]*1e8;
%Ht2 = [1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5 6.0 10 15 20]*1e8;
Hc = (Ht1+Ht2)/2;
Tc= 1e6*ones(2, length(Hc));
Tb= [10000; 10000];
Hb = 1.2e8;

param = reoGetParam;
param.wTemp = 3;
param.wTemp_D = 0;
param.rescntmax = 70;
reoIterationCore(hLib, mfoData, M, freqs, Robs, Lobs, H0, T0, Hb, Tb, Ht1, Ht2, Hc, Tc, NT, posR, posL, Mask, param, diagrH, diagrV, gridstep);

utilsFreeLibrary(hLib);

end


