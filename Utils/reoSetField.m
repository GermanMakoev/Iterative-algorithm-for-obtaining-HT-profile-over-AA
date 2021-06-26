function [MM, base, Bc] = reoSetField(hLib, mfoData, gridstep)

if ~exist('gridstep', 'var') || isempty(gridstep)
    gridstep = [1 1];
end

gridstepinR = gridstep/mfoData.R;
M = size(mfoData.B.z);
[MoutInt, base] = gstCalcSetField(hLib, ...
                  M, M, mfoData.B.x, mfoData.B.y, mfoData.B.z, mfoData.vcos, mfoData.stepP, mfoData.baseP, gridstepinR, mfoData.posangle);
base = base*mfoData.R; 
MM = double(MoutInt);

MM3 = [MM size(mfoData.B.z, 3)];
[Bc.x, Bc.y, Bc.z] = gstCalcGetField(hLib, MM3);

end
