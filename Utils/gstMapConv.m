function [scan, flux] = gstMapConv(map, diagr, diagrV, arcstep)
% diagrW in arcsec

sfumult = 2.35e8;

map(isnan(map)) = 0;
map = map .* repmat(diagrV', 1, size(map, 2));
Iint = sum(map, 1)*sfumult* arcstep(2);
scan = conv(Iint', diagr', 'same')';

flux = sum(scan)* arcstep(1);

end
