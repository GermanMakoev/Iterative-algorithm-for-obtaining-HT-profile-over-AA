function [diagrH, diagrV, diagrW, diagrWV] =  mfoCreateRATANDiagrams(freqs, sz, step, basearc, diagrMode)

lng = 3*sz(2);
diagrH = zeros([length(freqs) lng]);
diagrV = zeros([length(freqs) sz(1)]);
diagrW = zeros(1, length(freqs));
diagrWV = zeros(1, length(freqs));
%pts = linspace(-sz(1)/2*step(1), sz(1)/2*step(1), sz(1));
pts = basearc(1)/step(1) + (0:sz(1)-1);
for kf = 1:length(freqs) 
    diagrW(kf) = rtnGetDiagrammW(freqs(kf), diagrMode);
    diagrH(kf, :) = getGaussNorm(diagrW(kf)/step(2)/2, (lng+1)/2, lng);
    diagrWV(kf) = rtnGetDiagrammWV(freqs(kf));
    diagrV(kf, :) = getGaussPoints(diagrWV(kf)/step(1)/2, pts);
end

end
