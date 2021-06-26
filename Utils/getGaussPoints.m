function row = getGaussPoints(hw, points)

a0 = points(1);
M = (points(end) - a0)/(length(points)-1);
row = zeros(1, length(points));
m = -log(0.5)/hw^2;
for k = 1:size(row, 2)
    row(k) = exp(-(( (k-1)*M + a0 )^2)*m);
end

end
