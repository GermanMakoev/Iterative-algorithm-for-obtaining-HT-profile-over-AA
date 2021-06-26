function Mask = reoGetMaskA(B)

Mask = ones(size(B(:, :, 1)));
Mask(B(:, :, 1) > 1500) = 2;

end