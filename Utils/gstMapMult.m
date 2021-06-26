function map = gstMapMult(map0, diagr, diagrV, arcstep, pos)

map = map0;
map(isnan(map)) = 0;
map = map .* repmat(diagrV', 1, size(map, 2));
cd = floor((length(diagr)+1)/2);
xp = cd-pos+(1:size(map, 2));
dp = interp1(1:length(diagr), diagr, xp);
map = map .* repmat(dp, size(map, 1), 1);
map = map * 2.35e8*arcstep(1);

end
