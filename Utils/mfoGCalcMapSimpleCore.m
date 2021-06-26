function [pFOrd, pScanOrd, pFExt, pScanExt, base, rc] = mfoGCalcMapSimpleCore(filename, H, Temp, Dens, freq, harms, gridstep, maxB, vcos, posangle, freefree)

if ~exist('posangle', 'var') || isempty(posangle)
    posangle = 0;
end

if ~exist('gridstep', 'var') || isempty(gridstep)
    gridstep = [1 1];
end

if ~exist('freefree', 'var')
    freefree = 0;
end

[pFOrd, pScanOrd, pFExt, pScanExt] = deal([]);
rc = 0;

fullname = which(filename);
if isempty(fullname)
    disp('***** Field file does not exist *****');
    rc = 1;
    return
end

mfoData = iouLoadField(filename);
if isempty(mfoData) || ~isfield(mfoData, 'B') || isempty(mfoData.B)
    disp('***** Load field error *****');
    rc = 2;
    return
end

if exist('maxB', 'var') && ~isempty(maxB)
    B0 = sqrt(mfoData.B.x(:, :, 1).^2 + mfoData.B.y(:, :, 1).^2 + mfoData.B.z(:, :, 1).^2);
    mvB = max(max(abs(B0)));
    mfoData.B.x = mfoData.B.x/mvB*maxB;
    mfoData.B.y = mfoData.B.y/mvB*maxB;
    mfoData.B.z = mfoData.B.z/mvB*maxB;
end

if ~exist('vcos', 'var') || isempty(vcos)
    if isfield(mfoData, 'vcos')
        vcos = mfoData.vcos;
    else
        vcos = [0 0 1];
    end
end

libname = 'agsGeneralRadioEmission.dll';
fullname = which(libname);
if isempty(fullname)
    disp('***** Library file does not exist *****');
    return
end
hLib = utilsLoadLibrary(fullname);
if hLib == 0
    disp('***** Library load error *****');
    return
end
gstCalcInit(hLib);
gstSetPreferenceInt(hLib, 'cycloCalc.ConsiderFreeFree', freefree);

gridstepinR = gridstep/mfoData.R;
M = size(mfoData.B.z);
[MoutInt, base] = gstCalcSetField(hLib, ...
                  M, M, mfoData.B.x, mfoData.B.y, mfoData.B.z, vcos, mfoData.stepP, mfoData.baseP, gridstepinR, posangle);
base = base*mfoData.R; 
MM = double(MoutInt);
MM3 = [MM size(mfoData.B.z, 3)];

tic

mode = 2; %     enum Mode {USERDEF = 0, SIMPLE = 1, SINCALC = 2, APPROX = 3};
c = 0;
b = 0;
[pFOrd, pScanOrd, pFExt, pScanExt] = gstCalcMapScan(hLib, ...
                  length(H), H, Temp, Dens, ...
                  MM, freq, length(harms), harms, ...
                  mode, c, b, mfoData.R);

toc

mult = 2.35e8*gridstep(1)*gridstep(2);
pFOrd = pFOrd*mult;
pScanOrd = pScanOrd*mult;
pFExt = pFExt*mult;
pScanExt = pScanExt*mult;

utilsFreeLibrary(hLib);

end
