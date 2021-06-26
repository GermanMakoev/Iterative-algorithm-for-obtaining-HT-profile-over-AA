function [xsol, iter] = NelderMead(func, criteria, bound, x, reflex, expand, contract, shrink)

if ~exist('reflex', 'var')
    reflex = 1;
end
if ~exist('expand', 'var')
    expand = 2;
end
if ~exist('contract', 'var')
    contract = 0.5;
end
if ~exist('shrink', 'var')
    shrink = 0.5;
end

% initial simplex x(n+1, n), f(n+1, 1)
nsimp = size(x, 2) + 1;
f = zeros(nsimp, 1);
for k = 1:nsimp
    f(k) = feval(func, x(k, :));
end

iter = 0;
while ~feval(criteria, x, f)
    iter = iter + 1;
    fL = Inf;
    fG = 0;
    fH = 0;
    kH = [];
    kL = [];
    for k = 1:nsimp
        if f(k) > fH
            fH = f(k);
            kH = k;
        end
        if f(k) < fL
            fL = f(k);
            kL = k;
        end
    end
    for k = 1:nsimp
        if k == kH
            continue
        end
        if f(k) > fG
            fG = f(k);
        end
    end
    
    xL = x(kL, :);
    xH = x(kH, :);
    x0 = (sum(x, 1) - xH)/(nsimp-1);
    %f0 = feval(func, x0);
    
    xR = bound((1+reflex)*x0 - reflex*xH);
    fR = feval(func, xR);
    
    if fR < fL % good direction
        xE = bound(expand*xR + (1-expand)*x0);
        fE = feval(func, xE);
        if fE < fL % good expand, use it
            x(kH, :) = xE;
            f(kH) = fE;
        else  % bad expand, use reflexed
            x(kH, :) = xR;
            f(kH) = fR;
        end
    else % bad direction
        if fR < fG % reflexed is better than G, H
            x(kH, :) = xR;
            f(kH) = fR;
        else % collapse
            if fR < fH
                x(kH, :) = xR;
                f(kH) = fR;
            end
            xC = contract*xH + (1-contract)*x0;    
            fC = feval(func, xC);
            if fC < fH
                x(kH, :) = xC;
                f(kH) = fC;
            else % total collapse
                for k = 1:nsimp
                    x(k, :) = (x(k, :) + xL)*shrink;
                end
            end
        end
    end
end

xsol = xL;

end
