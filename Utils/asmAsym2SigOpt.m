function [appr, A, xc, w1, w2, w3] = asmAsym2SigOpt(freqs, fluxes)

smoo = fluxes;
smoo(2:end-1) = (smoo(1:end-2) + smoo(2:end-1) + smoo(3:end))/3;
[A, im] = max(smoo);
xc = freqs(im);
wingidxs = find(smoo >= A*0.9);
if isempty(wingidxs) || wingidxs(1) >= im
    frleft = freqs(im-2);
else
    frleft = (xc - freqs(1))/2;
end
if isempty(wingidxs) || wingidxs(end) <= im
    frright = freqs(im+2);
else
    frright = (xc + freqs(end))/2;
end
w1 = 10;
w2 = freqs(im) - frleft;
w3 = frright - freqs(im);

initsimp = [ ...
             A/2, xc-3, 0, w2, w3; ...
             A*2, xc-3, 0, w2, w3; ...
             A/2, xc+3, 0, w2, w3; ...
             A/2, xc-3, 5, w2, w3; ...
             A/2, xc-3, 0, w2+5, w3; ...
             A/2, xc-3, 0, w2, w3+5; ...
           ];


x = NelderMead(@l_calc, @l_crit, @l_bound, initsimp, 1, 2, 0.5, 0.5);

    function x = l_bound(x)
    end

    function f = l_calc(x)
        f = sum((fluxes - asmAsym2Sig(freqs, x(1), x(2), x(3), x(4), x(5))).^2);
    end

    function conv = l_crit(x, ~)
        conv = std(x(:, 1), 0, 1)/mean(x(:, 1), 1) < 1e-3 && all(std(x(:, 2:end), 0, 1) < 1e-3);
    end

A = x(1);
xc = x(2);
w1 = x(3);
w2 = x(4);
w3 = x(5);

appr = asmAsym2Sig(freqs, A, xc, w1, w2, w3);

end
