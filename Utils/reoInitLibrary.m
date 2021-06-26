function hLib = reoInitLibrary(QT, freefree, fullname)

if ~exist('fullname', 'var')
    libname = 'agsGeneralRadioEmission.dll';
    fullname = which(libname);
    if isempty(fullname)
        disp('***** Library file does not exist *****');
        return
    end
end

hLib = utilsLoadLibrary(fullname);
if hLib == 0
    disp('***** Library load error *****');
else
    gstCalcInit(hLib);
    
    if ~exist('QT', 'var')
        QT = 1;
    end
    gstSetPreferenceInt(hLib, 'cycloLine.ZoneSearch.QT.Use', QT);
    
    if ~exist('freefree', 'var')
        freefree = 0;
    end
    gstSetPreferenceInt(hLib, 'cycloCalc.ConsiderFreeFree', freefree);
end

end
