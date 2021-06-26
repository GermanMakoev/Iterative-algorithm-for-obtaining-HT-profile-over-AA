% загрузка магнитного поля (функция mfoLoadField заменила ранее используюмую iouLoadField)
mfoData = mfoLoadField('C:\Users\Alan Makoev\Desktop\Matlab — 2\SDO\12470_hmi.M_720s.20151218_082209.W85N13CR.CEA.NAS_1000_sst.sav'); % for az=+12, +10
if isempty(mfoData)
    return
end

% загрузка данных РАТАН (см. описание структуры файлов). Эта функция берет
% на себя чтение и загрузку данных в структуру ratan
ratan = iouLoadRATANData('C:\Users\Alan Makoev\Desktop\Matlab — 2\RATAN\RATAN_AR12470_20151218_082624_az+10_SPECTRA__flocculae-included_stille_appr.dat');

% в данном файле спектры в девяти точках скана, средняя точка (5) соответствует
% примерно максимуму поля, берем ее
select = 5;

% выделяем частоты и спектры
freqs = ratan.freqs;
Robs = ratan.right(select, :);
Lobs = ratan.left(select, :);

% можно работать и с исходными данными, но, чтобы убрать шум, можно
% попробовать аппроксимировать спектры подходящей функцией.
Robsappr = asmAsym2SigOpt(ratan.freqs*1e-9, Robs);
Lobsappr = asmAsym2SigOpt(ratan.freqs*1e-9, Lobs);
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

[M, base, Bc] = reoSetField(hLib, mfoData, step);

% данные РАТАН содержат позицию точек скана (в угл. сек.), поэтому позицию
% точки скана по горизонтали в терминах размеченного поля можно получить
% следующим образом:
pos = floor((ratan.pos(select) - base(2))/step(1)) + 1;

% там еще может быть сдвиг данных РАТАН, сейчас пока не включаем

% в данных SDO теперь имеются данные наблюдений в континууме, из которых
% можно пробовать строить маски для тени/полутени/флоккула, это тоже на
% следующем этапе.
