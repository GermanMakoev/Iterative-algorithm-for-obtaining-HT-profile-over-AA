function reoIterTest12470

mfoData = iouLoadField('C:\Users\Alan Makoev\Desktop\Matlab\12470_hmi.M_720s.20151218_085809.W84N12CR.CEA.NAS_sst_trim.mat'); % 50:240, 100:280
%mfoData = iouLoadField('s:\Projects\MLM\12470_hmi.M_720s.20151218_085809.W84N12CR.CEA.NAS_sst_trim.mat'); % 50:240, 100:280
if isempty(mfoData)
    return
end

QT = 1;
freefree = 0;

hLib = reoInitLibrary(QT, freefree);
if hLib == 0
    return
end

gridstep = [2 2];
[M, base, Bc] = reoSetField(hLib, mfoData, gridstep);
Mask = reoGetMaskA(sqrt(Bc.x.^2 + Bc.y.^2 + Bc.z.^2));

ratan = load('C:\Users\Alan Makoev\Desktop\Matlab\RATAN_12470_20151218_091001.mat'); %s:\University\Work\12470\ratan18_0.dat
%ratan = load('s:\Projects\MLM\RATAN_12470_20151218_091001.mat'); %s:\University\Work\12470\ratan18_0.dat
freqs=[ratan.freq];
[diagrH, diagrV] =  mfoCreateRATANDiagrams(freqs, M, gridstep, base, 3);
NT = 8e15;
H0 = [1, 1.4e8, 1.6e8, 5e8, 2e10];
T0 = [[1e4, 1e4, 2e6, 5e6, 5e6];[9e3, 1e4, 4e6, 5e6, 5e6]];
Dens= NT./T0;
%--------------------------------------------------------------------------
posR = 26*ones(1, length(freqs));
posL = 26*ones(1, length(freqs));
%-----------------------------------------------------------------------
Ht1 = [0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15 20 25]*1e8;
Ht2 = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 15 20 25 30]*1e8;
Hc = (Ht1+Ht2)/2;
Tc = 1e6*ones(2, length(Hc));
Tb = [6000 6000;6000 6000];
Hb = [0.7e8 0.8e8];

param = reoGetParam;
param.wL = 0.2;
param.wTemp = 30;
param.wTemp_D = 1e-7;
param.rescntmax = 30;
reoIterationCore1Point2Region(hLib, mfoData, M, freqs, [ratan.Robs], [ratan.Lobs], [], [], Hb, Tb, Ht1, Ht2, Hc, Tc, 3.2e15, posR, posL, Mask, param, diagrH, diagrV, gridstep);
utilsFreeLibrary(hLib);
end

