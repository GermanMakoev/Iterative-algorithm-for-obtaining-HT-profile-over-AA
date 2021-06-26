function map = ssuSymmColormap(nsteps, plus, minus, center, contrast)

if (~exist('nsteps', 'var'))
    nsteps = 231;
end;
if (~exist('plus', 'var'))
    plus = [1 0 0];
end
if (~exist('minus', 'var'))
    minus = [0 0 1];
end
if (~exist('center', 'var'))
    center = [0.9 0.9 0.9];
end
if (~exist('contrast', 'var'))
    contrast = 3;
end

zero = [1 1 1];

nsteps = floor(nsteps/2)*2 + 1;
n2 = floor(nsteps/2)+1;
mult = ((0:(n2-1))/(n2-1)).^contrast;
multmap = repmat(mult', 1, 3);

map = zeros(nsteps, 3);

dminus = (center - minus);
minusmap = repmat(dminus, n2, 1);
map(1:n2, :) = minusmap .* multmap + repmat(minus, n2, 1);

dplus = (center - plus);
plusmap = repmat(dplus, n2, 1);
map(n2:end, :) = plusmap .* flipdim(multmap, 1) + repmat(plus, n2, 1);

map(n2, :) = zero;

end
