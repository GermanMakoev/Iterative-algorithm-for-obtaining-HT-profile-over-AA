function mfoData = iouLoadField(filename)

mfoData = [];

if ~exist('filename', 'var')
    filename = wfr_getFieldFilename;
end

if ~exist(filename, 'file')
    filename = which(filename);
end

if isempty(filename)
    disp(['***** Field file "' filename '" does not exist *****']);
    return
end

mfoData.filename = filename;
mfoData.filetype = '';
mfoData.fileid = '';
mfoData.type = '';
mfoData.B = [];
mfoData.Bpot = [];
mfoData.Bmod = [];
mfoData.time = '';
mfoData.lon = [];       % degrees
mfoData.lat = [];       % degrees
mfoData.step = [];      % arcsec x 3
mfoData.centre = [];    % arcsec x 2
mfoData.stepP = [];     % R=1 x 3
mfoData.baseP = [];     % R=1 x 2
mfoData.R = [];         % arcsec
mfoData.vcos = [];

try
    s = load(filename, '-mat');
    if isfield(s, 'trim') && isfield(s.trim, 'Result') && isfield(s.trim.Result, 'LFFF') % HMI
        mfoData.filetype = 'HMI';
        B = fpuField2XYZ(s.trim.Result.LFFF.B4);
        if isfield(s.trim.Result, 'Wiegelmann')
            mfoData.type = 'NLFFF';
            mfoData.B = fpuField2XYZ(s.trim.Result.Wiegelmann.B4);
            mfoData.Bpot = B;
        else
            mfoData.type = 'LFFF';
            mfoData.B = B;
        end
        if isfield(s.data, 'model')
            mfoData.Bmod = s.data.model.B;
        end
        mfoData.time = [s.OBS_DATE s.OBS_TIME];
        mfoData.step = [s.trim.Result.transform.arcstep(1) s.trim.Result.transform.arcstep(2) s.trim.Result.transform.arcstep(1)];
        mfoData.centre = [s.trim.YCEN s.trim.XCEN];
        mfoData.R = s.SOLAR_R;         % arcsec
        cp = mfoData.centre/mfoData.R;
        mfoData.lat = asind(cp(1));
        coslat = sqrt(1-cp(1).*cp(1));
        mfoData.lon = asind(cp(2)/coslat);
        mfoData.stepP = mfoData.step/mfoData.R;
        mfoData.baseP = l_base(mfoData.B.x, mfoData.stepP(1));
        mfoData.vcos = l_vcos(mfoData.lon, mfoData.lat);
    elseif isfield(s, 'B') && isfield(s, 'modstep')
        mfoData.filetype = 'MAT';
        mfoData.type = 'MOD';
        mfoData.B = s.B;
        mfoData.stepP = [s.modstep s.modstep s.modstep];
        mfoData.baseP = l_base(mfoData.B.x, mfoData.stepP(1));
        mfoData.R = 950;
        mfoData.step = mfoData.stepP*mfoData.R;
        mfoData.vcos = [0 0 1];
    elseif isfield(s, 'mfoData')
        mfoData = s.mfoData;
    else
        disp(['***** Wrong MAT format "' filename '" *****'])
        mfoData = [];
    end
catch sle % not mat-file, try others
    out = iouSAV2Data(filename);
    if isfield(out, 'MFODATA')
        out = out.MFODATA;
    end
    if ~isempty(out) && isfield(out, 'BZ') % SAV
        mfoData.filetype = 'SAV';
        if strfind(out.FILEID, '.POT')
            mfoData.type = 'LFFF';
        else
            mfoData.type = 'NLFFF';
        end
        mfoData.B.x = out.BX;
        mfoData.B.y = out.BY;
        mfoData.B.z = out.BZ;
        if isfield(out, 'BZPOT')
            mfoData.Bpot.x = out.BXPOT;
            mfoData.Bpot.y = out.BYPOT;
            mfoData.Bpot.z = out.BZPOT;
        end
        if isfield(out, 'TIME')
            mfoData.time = out.TIME;
        elseif isfield(out, 'OBSTIME')
            mfoData.time = out.OBSTIME;
        end
        mfoData.fileid = out.FILEID;
        mfoData.lon = out.LONV;
        mfoData.lat = out.LATV;
        mfoData.step = [out.D_ARC out.D_ARC out.D_ARC];
        mfoData.centre = out.CENTRE';
        mfoData.R = out.R_ARC;
        
        mfoData.stepP = mfoData.step/mfoData.R;
        mfoData.baseP = l_base(mfoData.B.x, mfoData.stepP(1));
        
        mfoData.vcos = l_vcos(mfoData.lon, mfoData.lat);
    else
        out = iouBIN2Data(filename);
        if ~isempty(out) && isfield(out, 'BZ') % BIN/OUT
            mfoData.filetype = 'BIN';
            if ~isfield(out, 'XBOUNDS')
                mfoData.type = 'LFFF';
            else
                mfoData.type = 'NLFFF';
            end
            mfoData.B.x = out.BX;
            mfoData.B.y = out.BY;
            mfoData.B.z = out.BZ;
        else
            disp(['***** Wrong non-MAT format "' filename '" *****'])
            mfoData = [];
        end
    end
end



end

%--------------------------------------------------------------------------
function base = l_base(b, modstep)

[nx, ny, ~] = size(b);
base = [-modstep*0.5*(nx-1), -modstep*0.5*(ny-1)];

end

%--------------------------------------------------------------------------
function vcos = l_vcos(lon, lat)

vcos(1) = -cosd(lon)*sind(lat);
vcos(2) = -sind(lon);
vcos(3) =  cosd(lon)*cosd(lat);

end

