function row = getGaussNorm(hw, pos, lng)

row = zeros(1, lng);
m = -log(0.5)/hw^2;
for k = 1:lng
    row(k) = exp(-((k-pos)^2)*m);
end

row = sqrt(m/pi)*row;

end
