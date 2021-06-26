function data = iouLoadRATANData(filename)

[header, params, row, mode] = l_readHeader(filename);
M = dlmread(filename, '', row, 0);

if strcmp(mode, 'spectra')
    data = struct('mode', mode, 'header', header, 'params', params, 'pos', [], 'freqs', [], 'right', [], 'left', [], 'qs', [], 'ffpos', [], 'ff', []);
    data.header.shift = M(1, 1);
    nflux = data.header.N_POS * 3;
    data.pos = (M(1, (1:3:nflux)+1))';
    data.freqs = (M(2:end, 1)*1e9)';
    if nflux+1 < size(M, 2)
        data.ffpos = M(1, nflux+2);
        data.ff    = (M(2:end, nflux+2)*1e-4)';
    end
else
    data = struct('mode', mode, 'header', header, 'params', params, 'pos', [], 'freqs', [], 'right', [], 'left', [], 'qs', []);
    data.header.shift = M(1, 1);
    nflux = data.header.N_FREQS * 3;
    data.freqs = M(1, (1:3:nflux)+1)*1e9;
    data.pos   = M(2:end, 1);
end
data.right = M(2:end, (1:3:nflux)+1)*1e-4;
data.left  = M(2:end, (1:3:nflux)+2)*1e-4;
data.qs    = M(2:end, (1:3:nflux)+3)*1e-4;

if strcmp(mode, 'spectra')
    data.right = data.right';
    data.left = data.left';
    data.qs = data.qs';
end

end

%--------------------------------------------------------------------------
function [header, params, row, mode] = l_readHeader(filename)

row = 0;
header = [];
params = [];

mode = '';
fid = fopen(filename, 'r');

predefined = {'RATAN Spectra-At-Positions Data File ', 'version';
              'RATAN Selected Scans Data File ', 'version';
              'generated at ', 'created';
              'frequencies selected by user for scan', 'user_sel_freqs'
              };
while (~feof(fid))
    sline = fgets(fid);

    if ~strcmp(sline(1), '#')
        break
    end
    sline = strtrim(sline(3:end));
    row = row + 1;
    
    predef = false;
    for k = 1:size(predefined, 1)
        lng = length(predefined{k, 1});
        if length(sline) > lng && strcmp(sline(1:lng), predefined{k, 1})
            header.(predefined{k, 2}) = strtrim(sline((lng+1):end-1));
            predef = true;
            if k == 1
                mode = 'spectra';
            elseif k == 2
                mode = 'scans';
            end
            
            break
        end
    end
    
    if predef
        continue
    end
    
    outs = regexp(sline, '([^ ]*) += +(.*)', 'tokens');
    name = outs{1}{1};
    if any(strcmp(outs{1}{1}, {'DATE-OBS', 'TIME-OBS'}))
        v = strtrim(outs{1}{2});
    else
        v = l_getval(outs{1}{2});
    end
    name = strrep(name, '-', '_');
    par = regexp(name, 'PAR\.(.*)', 'tokens');
    if isempty(par)
        header.(name) = v;
    else
        params.(par{1}{1}) = v;
    end
end

fclose(fid);

end

%--------------------------------------------------------------------------
function v = l_getval(s)

v = str2num(s);
if isempty(v)
    v = strtrim(s);
end

end

