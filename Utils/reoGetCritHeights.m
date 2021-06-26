function [heights, Bcrits, Bmax] = reoGetCritHeights(mfoData, freqs, harm)

B = fpuFieldVal(mfoData.B);
Bmax = squeeze(max(max(B, [], 1), [], 2));
HB = mfoData.step(3)/mfoData.R*6.96e10*(0:size(mfoData.B.z, 3)-1);

heights = zeros(1, length(freqs));
Bcrits = zeros(1, length(freqs));
for k = 1:length(freqs)
    Bcrits(k) = freqs(k)/harm/2.799e6;
    heights(k) = interp1(Bmax, HB, Bcrits(k));
end

end
