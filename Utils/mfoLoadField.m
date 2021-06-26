function [mfoData, filename, out] = mfoLoadField(filename) %, rangev, rangeh)

mfoData = [];
out = [];

if ~exist('filename', 'var') || length(filename) == 3
    if exist('filename', 'var') && length(filename) == 3
        ext = filename;
    else
        ext = '';
    end
    filename = wfr_getFieldFilename(ext);
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
mfoData.version = '';
mfoData.version_info = '';
mfoData.B = [];
mfoData.Bpot = [];
mfoData.Bmod = [];
mfoData.cont = [];
mfoData.mask = [];
mfoData.time = '';
mfoData.lon = 0;       % degrees
mfoData.lat = 0;       % degrees
mfoData.step = [1 1 1];      % arcsec x 3
mfoData.centre = [0 0];    % arcsec x 2
mfoData.stepP = [];     % R=1 x 3
mfoData.baseP = [];     % R=1 x 2
mfoData.R = 960;         % arcsec
mfoData.vcos = [0 0 1];
mfoData.AIA = [];

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
            mfoData.Bpot = B;
        end
        if isfield(s.data, 'model') && isfield(s.data.model, 'B')
            mfoData.Bmod = s.data.model.B;
        end
        mfoData.time = [s.OBS_DATE s.OBS_TIME];
        step = s.trim.Result.transform.arcstep(1);
        step2 = step;
        if length(s.trim.Result.transform.arcstep) > 1
            step2 = s.trim.Result.transform.arcstep(2);
        end
        mfoData.step = [step step2 step];
        mfoData.centre = [s.trim.YCEN s.trim.XCEN];
        mfoData.R = s.SOLAR_R;         % arcsec
        [mfoData.lat, mfoData.lon] = l_latlon(mfoData.centre/mfoData.R);
        mfoData.stepP = mfoData.step/mfoData.R;
        mfoData.baseP = l_base(mfoData);
        mfoData.vcos = crdVCos(mfoData.lat, mfoData.lon);
        mfoData.trim = s.trim;
    elseif isfield(s, 'B') && isfield(s, 'modstep')
        mfoData.filetype = 'MAT';
        mfoData.type = 'MOD';
        mfoData.B = s.B;
        mfoData.stepP = [s.modstep s.modstep s.modstep];
        mfoData.baseP = l_base(mfoData);
        mfoData.step = mfoData.stepP*mfoData.R;
    elseif isfield(s, 'mfoData')
        mfoData.filetype = 'MAT';
        mfoData = s.mfoData;
    elseif isfield(s, 'BZ')
        mfoData.filetype = 'PRP';
        mfoData.B.x = s.BX;
        mfoData.B.y = s.BY;
        mfoData.B.z = s.BZ;
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
        if ~isfield(out, 'SST_VERSION')
            mfoData = l_mfosavprev(mfoData, out);
        elseif strcmp(out.SST_VERSION, '20200228')
            mfoData = l_mfosav20200228(mfoData, out);
        else
            disp(['***** sst format "' out.SST_VERSION '" not implemented. *****'])
            mfoData = [];
        end    
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

% if ~isempty(mfoData) && (exist('rangev', 'var') || exist('rangeh', 'var'))
%     if ~exist('rangev', 'var') || isempty(rangev)
%         rangev = 1:size(mfoData.B.z, 1);
%     end
%     if ~exist('rangeh', 'var') || isempty(rangeh)
%         rangeh = 1:size(mfoData.B.z, 2);
%     end
%     
%     mfoData = mfoTrimData(mfoData, rangev, rangeh);
% end

end

%--------------------------------------------------------------------------
function base = l_base(mfoData)

[nx, ny, ~] = size(mfoData.B.x);
%base = mfoData.centre/mfoData.R - [mfoData.stepP(1)*0.5*(nx-1), mfoData.stepP(2)*0.5*(ny-1)];
base =  - [mfoData.stepP(1)*0.5*(nx-1), mfoData.stepP(2)*0.5*(ny-1)];

end

%--------------------------------------------------------------------------
function [lat, lon] = l_latlon(cp)

lat = asind(cp(1));
coslat = sqrt(1-cp(1).*cp(1));
lon = asind(cp(2)/coslat);

end

%--------------------------------------------------------------------------
function mfoData = l_mfosavcommon(mfoData, out)

if ~isfield(out, 'FILEID')
    mfoData.type = 'None';
else
    if strfind(out.FILEID, '.POT')
        mfoData.type = 'LFFF';
    elseif strfind(out.FILEID, '.BND')
        mfoData.type = 'BOUNDED';
    elseif strfind(out.FILEID, '.NAS')
        mfoData.type = 'NLFFF';
    else
        mfoData.type = 'Unknown';
    end
end

mfoData.B.x = out.BX;
mfoData.B.y = out.BY;
mfoData.B.z = out.BZ;
mfoData.cont = getfieldsafe(out, 'IC', []);
if isfield(out, 'TIME')
    mfoData.time = out.TIME;
elseif isfield(out, 'OBSTIME')
    mfoData.time = out.OBSTIME;
end
mfoData.fileid = getfieldsafe(out, 'FILEID', '');

mfoData.R = out.R_ARC;                              % arcsec

mfoData.stepP = mfoData.step/mfoData.R;
mfoData.baseP = l_base(mfoData);
if isfield(out, 'AIA_IDS') && ~isempty(out.AIA_IDS)
    mfoData.AIA.Ids = out.AIA_IDS;
    mfoData.AIA.size = out.AIA_SIZE;
    mfoData.AIA.center = out.AIA_CENTER;
    mfoData.AIA.step = out.AIA_STEP;
    mfoData.AIA.RSun = out.AIA_RSUN;
    mfoData.AIA.data = out.AIA_DATA;
end
if isfield(out, 'MODEL_MASK') && ~isempty(out.MODEL_MASK)
    mfoData.mask = double(out.MODEL_MASK);
end
end

%--------------------------------------------------------------------------
function mfoData = l_mfosavprev(mfoData, out)

if isfield(out, 'BZPOT')
    mfoData.Bpot.x = out.BXPOT;
    mfoData.Bpot.y = out.BYPOT;
    mfoData.Bpot.z = out.BZPOT;
end
mfoData.lon = getfieldsafe(out, 'LONV', 0);
mfoData.lat = getfieldsafe(out, 'LATV', 0);
if isfield(out, 'D_ARC')
    mfoData.step = [out.D_ARC out.D_ARC out.D_ARC];
end
mfoData.centre = getfieldsafe(out, 'CENTRE', mfoData.centre);
mfoData.centre = getfieldsafe(out, 'CENTER', mfoData.centre);
mfoData.centre = mfoData.centre';
mfoData.R = getfieldsafe(out, 'R_ARC', mfoData.R);
mfoData.lon_vis = getfieldsafe(out, 'LON_VIS', mfoData.R);
mfoData.lat_vis = getfieldsafe(out, 'LAT_VIS', mfoData.R);

mfoData.vcos = crdVCos(mfoData.lat, mfoData.lon);

mfoData = l_mfosavcommon(mfoData, out);

end

%--------------------------------------------------------------------------
function mfoData = l_mfosav20200228(mfoData, out)

mfoData.version = '20200228';
mfoData.version_info = out.VERSION_INFO;

mfoData.lon = out.LON_CEN;       % degrees
mfoData.lat = out.LAT_CEN;       % degrees

mfoData.step = [out.DX_ARC out.DY_ARC out.DX_ARC];  % arcsec x 3
mfoData.centre = [out.X_CEN out.Y_CEN]*out.R_ARC;   % arcsec x 2
mfoData.x_box = out.X_BOX;
mfoData.y_box = out.Y_BOX;
mfoData.lon_hg = out.LON_HG;
mfoData.lat_hg = out.LAT_HG;

mfoData.vcos = out.VCOS;
if isfield(out, 'AIA') && isfield(out.AIA, 'IDS') && ~isempty(out.AIA.IDS)
    mfoData.AIA.Ids = out.AIA_IDS;
    mfoData.AIA.size = out.AIA_SIZE;
    mfoData.AIA.center = out.AIA_CENTER;
    mfoData.AIA.step = out.AIA_STEP;
    mfoData.AIA.RSun = out.AIA_RSUN;
    mfoData.AIA.data = out.AIA_DATA;
end

mfoData = l_mfosavcommon(mfoData, out);

end

