function reoIterTest12794v4

% загрузка магнитного поля (функция mfoLoadField заменила ранее используемую iouLoadField)
mfoData = mfoLoadField('SDO\12794_hmi.M_720s.20201228_091030.E16S16CR.CEA.NAS_750_sst.sav'); % for az=+12, +10
if isempty(mfoData)
    return
end

% загрузка данных РАТАН (см. описание структуры файлов). Эта функция берет
% на себя чтение и загрузку данных в структуру ratan
ratan = iouLoadRATANData('RATAN\RATAN_AR12794_20201228_091520_az0_SPECTRA.dat');

% в данном файле спектры в девяти точках скана, средняя точка (5) соответствует
% примерно максимуму поля, берем ее
select = [2,3,4];

% выделяем частоты и спектры
freqs = ratan.freqs;
Robs = zeros(length(freqs),3);
Lobs = zeros(length(freqs),3);
Robsappr = zeros(length(freqs),3);
Lobsappr = zeros(length(freqs),3);
posR = zeros(length(freqs),3);
posL = zeros(length(freqs),3);
Points = length(select); 
for l=1:3
%Robs(:,l) = ratan.left(select(l), :);
%Lobs(:,l) = ratan.right(select(l), :);
Robs(:,l) = ratan.left(select(l), 1:length(freqs));
Lobs(:,l) = ratan.right(select(l), 1:length(freqs));
end 
% можно работать и с исходными данными, но, чтобы убрать шум, можно
% попробовать аппроксимировать спектры подходящей функцией.
for l=1:3
Robsappr(:,l) = asmAsym2SigOpt(ratan.freqs*1e-9, Robs(:,l)');
Lobsappr(:,l) = asmAsym2SigOpt(ratan.freqs*1e-9, Lobs(:,l)');
% для уверенности в том, что аппроксимация приемлема, пожно построить,
% скажем, Robs и Robsappr на одном графике
end 

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

%Увеличение поля
%FF=-0.2;
%mfoData.B.x=mfoData.B.x*(1+FF);
%mfoData.B.y=mfoData.B.y*(1+FF);
%mfoData.B.z=mfoData.B.z*(1+FF);

step=[2 2];
[M, base, Bc] = reoSetField(hLib, mfoData, step);
%Mask = mfoData.mask-5;
load('Mask12794.mat','Mask');
[diagrH, diagrV] =  mfoCreateRATANDiagrams(freqs, M, step, base, 3);

% данные РАТАН содержат позицию точек скана (в угл. сек.), поэтому позицию
% точки скана по горизонтали в терминах размеченного поля можно получить
% следующим образом:
for l=1:3
posR(:,l) = floor((ratan.pos(select(l)) - base(2))/step(1)) + 1;
posL(:,l) = floor((ratan.pos(select(l)) - base(2))/step(1)) + 1;
Points(l) = floor((ratan.pos(select(l)) - base(2))/step(1)) + 1;
end
% там еще может быть сдвиг данных РАТАН, сейчас пока не включаем

% в данных SDO теперь имеются данные наблюдений в континууме, из которых
% можно пробовать строить маски для тени/полутени/флоккула, это тоже на
% следующем этапе.

%-----------------------------------------------------------------------
Ht1 = [0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15 20 25]*1e8;
Ht2 = [1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 15 20 25 30]*1e8;
%Ht1 = [1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 8 10  12 15 20 25]*1e8;
%Ht2 = [1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 12 15 17 20 25 30]*1e8;
Hc = (Ht1+Ht2)/2;
Tc = 1e6*ones(2, length(Hc));
%Tc(1,1:10) = 6000;
Tb = [6000, 6000; 6000, 6000];

%Hb = [1.2e8, 1.3e8];
Hb = [0.8e8, 0.9e8];
%NT = [1e15, 2e15, 3e15, 4e15, 5e15, 1e16, 2e16, 5e16];
NT = 1e16;
%NT = 5e16;

param = reoGetParam;
param.wTemp = 1000;
%param.wL = 0.2;
%param.wL = 0.005;
param.rescntmax = 30;
reoIterationCore3Points2RegionsV2(hLib, mfoData, M, freqs, Robsappr, Lobsappr, [], [], Hb, Tb, Ht1, Ht2, Hc, Tc, NT, posR, posL, Mask, param, diagrH, diagrV, step, Points);
utilsFreeLibrary(hLib);
end

