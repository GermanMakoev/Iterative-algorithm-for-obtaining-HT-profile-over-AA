function reoIterTest12419v2

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
select = [2,5,6];

% выделяем частоты и спектры
freqs = ratan.freqs;
Robs = zeros(length(freqs),3);
Lobs = zeros(length(freqs),3);
Robsappr = zeros(length(freqs),3);
Lobsappr = zeros(length(freqs),3);
posR = zeros(length(freqs),3);
posL = zeros(length(freqs),3);
Points = zeros(length(select));
for l=1:3
Robs(:,l) = ratan.right(select(l), :);
Lobs(:,l) = ratan.left(select(l), :);
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

step=[2,2];
[M, base, Bc] = reoSetField(hLib, mfoData, step);

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

freqs(58:end) = [];
Robs(58:end,:) = [];
Lobs(58:end,:) = [];

%-----------------------------------------------------------------------
Ht1 = [0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.5 4.0 4.5 5.0 5.5  6 10 15 20 25]*1e8;
Ht2 = [1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.3 2.5 2.7 2.9 3.1 3.7 4.2 4.7 5.2 5.7 6.5 10 15 20 25 30]*1e8;
Hc = (Ht1+Ht2)/2;
Tc = 1e6*ones(1, length(Hc));
Tb = 6000;
Hb = 0.8e8;
NT = 2e15;

param = reoGetParam;
param.wTemp = 300;
param.wL = 0.2;
param.rescntmax = 30;
reoIterationCore3Points1Region(hLib, mfoData, M, freqs, Robs, Lobs, [], [], Hb, Tb, Ht1, Ht2, Hc, Tc, NT, posR, posL, param, Points);
utilsFreeLibrary(hLib);
end

