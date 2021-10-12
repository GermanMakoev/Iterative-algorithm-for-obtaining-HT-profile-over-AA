function reoIterTest12419v3

% загрузка магнитного поля (функция mfoLoadField заменила ранее используюмую iouLoadField)
mfoData = mfoLoadField('SDO\12419_hmi.M_720s.20150918_094619.E164N13CR.CEA.NAS_400_sst.sav'); % for az=+12, +10
if isempty(mfoData)
    return
end

% загрузка данных РАТАН (см. описание структуры файлов). Эта функция берет
% на себя чтение и загрузку данных в структуру ratan
ratan = iouLoadRATANData('RATAN\RATAN_AR12419_20150918_095555_az-12_SPECTRA__lin_appr.dat');

% в данном файле спектры в девяти точках скана, средняя точка (5) соответствует
% примерно максимуму поля, берем ее
select = 6;

% выделяем частоты и спектры
freqs = ratan.freqs;
Robs = ratan.right(select, :);
Lobs = ratan.left(select, :);

% можно работать и с исходными данными, но, чтобы убрать шум, можно
% попробовать аппроксимировать спектры подходящей функцией.
%Robsappr = asmAsym2SigOpt(ratan.freqs*1e-9, Robs);
%Lobsappr = asmAsym2SigOpt(ratan.freqs*1e-9, Lobs);
% для уверенности в том, что аппроксимация приемлема, пожно построить,
% скажем, Robs и Robsappr на одном графике

% позиционный угол РАТАН уже содержится в данных, поэтому просто
% перетаскиваем его в mfoData
mfoData.posangle = ratan.header.RATAN_P;

% далее (пропущено) обычная загрузка и инициализация библиотеки, в частности, разметка
% поля:
QT = 1;
freefree = 0;

hLib = reoInitLibrary(QT, freefree);
if hLib == 0
    return
end

step=[2 2];
[M, base, Bc] = reoSetField(hLib, mfoData, step);
%Mask = mfoData.mask;
load('Mask12419.mat','Mask');
[diagrH, diagrV] =  mfoCreateRATANDiagrams(freqs, M, step, base, 3);

% данные РАТАН содержат позицию точек скана (в угл. сек.), поэтому позицию
% точки скана по горизонтали в терминах размеченного поля можно получить
% следующим образом:
posR(1:length(freqs)) = floor((ratan.pos(select) - base(2))/step(1)) + 1;
posL(1:length(freqs)) = floor((ratan.pos(select) - base(2))/step(1)) + 1;
%posR(1:length(freqs))=48;
%posL(1:length(freqs))=48;
% там еще может быть сдвиг данных РАТАН, сейчас пока не включаем

% в данных SDO теперь имеются данные наблюдений в континууме, из которых
% можно пробовать строить маски для тени/полутени/флоккула, это тоже на
% следующем этапе.

freqsR = freqs(1:25);
freqsL = freqs(1:57);
Robs(26:end) = [];
Lobs(58:end) = [];

%-----------------------------------------------------------------------
Ht1 = [0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15 20 25]*1e8;
Ht2 = [1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 15 20 25 30]*1e8;
Hc = (Ht1+Ht2)/2;
Tc = 1e6*ones(2, length(Hc));
%Tc(2,1:4) = 6000;
Tb = [6000, 6000; 6000, 6000];
Hb = [0.8e8, 0.9e8];
NT = 2e16;

param = reoGetParam;
param.wTemp = 300;
param.wL = 0.2; 
param.rescntmax = 30;
reoIterationCore1Point2RegionV2(hLib, mfoData, M, freqs, freqsR, freqsL, Robs, Lobs, [], [], Hb, Tb, Ht1, Ht2, Hc, Tc, NT, posR, posL, Mask, param, diagrH, diagrV, step);
utilsFreeLibrary(hLib);
end

