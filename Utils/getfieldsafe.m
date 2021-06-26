function [v, found] = getfieldsafe(st, fld, dflt)

found = false;

% Returns requested field of a struct, or the default value, if the field is missing
if isfield(st, fld)
    v = [st.(fld)];
    found = true;

elseif nargin > 2
    v = dflt;

else
    v = [];
end

end