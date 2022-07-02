function reoIterationCore1of3Point1Region(hLib, mfoData, Mmap, freqs, Robs, Lobs, RB, LB, H0, T0, Hb, Tb, Ht1, Ht2, Hc, Tc, NTv, posR, posL, param)

nNT = 1;
NT = NTv(nNT);

outCntres = zeros(length(NTv));
outCntNT = zeros(length(NTv));
outRes = zeros(length(NTv));

[~, hAxes] = ssuCreateMultiAxesFig(2, 2, '', false, false);
hT0 = [];
title(hAxes(1,1),'ТВП');
xlabel(hAxes(1,1),'Высота, см');
ylabel(hAxes(1,1),'Температура, К');
title(hAxes(1,2),'Зависимость невязки от числа итераций');
xlabel(hAxes(1,2),'Число итераций');
ylabel(hAxes(1,2),'Невязка');
title(hAxes(2,1),'Спектры, правая поляризация');
xlabel(hAxes(2,1),'Частота, Гц');
ylabel(hAxes(2,1),'Поток, с.е.п./угл.сек.');
title(hAxes(2,2),'Спектры, левая поляризация');
xlabel(hAxes(2,2),'Частота, Гц');
ylabel(hAxes(2,2),'Поток, с.е.п./угл.сек.');
if ~isempty(H0) && ~isempty(T0)
    hT0 = plot(hAxes(1,1), H0, T0, 'Color', 'g');
end
set(hAxes(1,1), 'YScale', 'log', 'XLim', [0 1e9])
set(hAxes(1,2), 'YScale', 'log')
set(hAxes,'Fontsize',16)
hTc = [];
hRs0 = plot(hAxes(2,1), freqs, Robs, 'Color', 'g');
hRsc = [];
hLs0 = plot(hAxes(2,2), freqs, Lobs, 'Color', 'g');
hLsc = [];
drawnow

resRv = [];
resLv = [];
Tprev = Tc;
relax = param.relax;
cntres = 0;
cntNT = 0;
resprev = 0;
x = ones(1, length(Hc));
while true    
    cntR = zeros(length(freqs),length(Hc));
    cntL = zeros(length(freqs),length(Hc));
    
    Tc = relax*x.*Tc + (1-relax)*Tprev;
    Tc(Tc < param.Tmin) = param.Tmin;
    Tprev = Tc;
    
    Hcalc = [Hb Hc];
    Tcalc = [Tb Tc];
    Dcalc = NT./Tcalc;
    
    R = zeros(length(freqs), length(Hc));
    L = zeros(length(freqs), length(Hc));
    HRm = zeros(length(freqs), length(Hc));
    HLm = zeros(length(freqs), length(Hc));
    
    Rcalc = zeros(1, length(freqs));
    Lcalc = zeros(1, length(freqs));
    for k = 1:length(freqs)      
            posRcalc = posR(k);
            posLcalc = posL(k);
        
        [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
         depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft, ...
         pScanRight, pScanLeft] = ...
                 reoCalculate(hLib, mfoData, Hcalc, Tcalc, Dcalc, Mmap, freqs(k), param.harms, param.pTauL, param.mode, param.c, param.b, posRcalc, posLcalc);
        Rcalc(k) = pScanRight(posR(k));
        Lcalc(k) = pScanLeft(posL(k));
                 
        for i= 1:size(pFRight, 1)
            for j = 1:size(pFRight, 2)
                for m = 1:length(Hc)
                    if depthRight(i,j) >= 1 && Ht1(m) <= pHLRight(i,j,1) && pHLRight(i,j,1) < Ht2(m)
                        R(k,m) = R(k,m) + pFLRight(i,j,1);
                        cntR(k,m) = cntR(k,m) + 1;
                        HRm(k,m) = HRm(k,m) + pHLRight(i,j,1);
                    end
                    if depthLeft(i,j) >= 1 && Ht1(m) <= pHLLeft(i,j,1) && pHLLeft(i,j,1) < Ht2(m)
                        L(k,m) = L(k,m) + pFLLeft(i,j,1);
                        cntL(k,m) = cntL(k,m) + 1;
                        HLm(k,m) = HLm(k,m) + pHLLeft(i,j,1);
                    end
                end
            end
        end
        HRm(k,:) = HRm(k,:)./cntR(k,:);
        HLm(k,:) = HLm(k,:)./cntL(k,:);
    end
    
Rcalc = Rcalc + RB;
Lcalc = Lcalc + LB;

    resR = sqrt(sum((Robs-Rcalc).^2));
    resRrel = resR/sum(Robs);
    resL = sqrt(sum((Lobs-Lcalc).^2));
    resLrel = resL/sum(Lobs);
    resRv = [resRv resRrel];
    resLv = [resLv resLrel];
    resS = resRrel + resLrel;
    if abs(resS - resprev)/resS < param.reslim
        cntres = cntres + 1;
    else
        cntres = 0;
    end
    if cntres >= param.rescntstab || cntNT >= param.rescntmax
        disp(['*NT=' num2str(NT) ' cntStable=' num2str(cntres) ' cntNT=' num2str(cntNT) ' res=' num2str(resS)]);
        outCntres(nNT) = cntres;
        outCntNT(nNT) = cntNT;
        outRes(nNT) = resS;
        
        if nNT >= length(NTv)
            break
        end
        nNT = nNT + 1;
        NT = NTv(nNT);
        cntres = 0;
        cntNT = 0;
    end
    resprev = resS;
    cntNT = cntNT + 1;
    disp(['  cntStable=' num2str(cntres) ' cntNT=' num2str(cntNT) ' res=' num2str(resS)]);
    
    [hT0, hTc] = l_updatePlot(hAxes(1, 1), hT0, hTc, H0, T0, Hcalc, Tcalc, 'm');
    [hRs0, hRsc] = l_updatePlot(hAxes(2, 1), hRs0, hRsc, freqs, Robs, freqs, Rcalc, 'r');
    [hLs0, hLsc] = l_updatePlot(hAxes(2, 2), hLs0, hLsc, freqs, Lobs, freqs, Lcalc, 'b');
    lr = length(resRv);
    if lr > 1
        plot(hAxes(1,2), [lr-1 lr], resRv(end-1:end), 'Color', 'r')
        plot(hAxes(1,2), [lr-1 lr], resLv(end-1:end), 'Color', 'b')
        plot(hAxes(1,2), [lr-1 lr], resRv(end-1:end)+resLv(end-1:end), 'Color', 'k')
    end
    drawnow
    
    inclR = true(1, length(freqs));
    for i = 1:length(freqs)
        if all(R(i, :) == 0)
            inclR(i) = false;
        end
    end
    inclL = true(1, length(freqs));
    for i = 1:length(freqs)
        if all(L(i, :) == 0)
            inclL(i) = false;
        end
    end
    
    M = [R(inclR, :); L(inclL, :)];
    Fobs =  [Robs(inclR)  Lobs(inclL)];
    Fcalc = [Rcalc(inclR) Lcalc(inclL)];
    Fr = zeros(1, length(Fobs));
    for i = 1:length(Fobs)
        Fr(i) = Fobs(i) + sum(M(i, :)) - Fcalc(i);
    end
    sm = sum(sum(M))/numel(M(M > 0))*1e-6;
    nFreqEq = length(Fr);
    
    Z = zeros(length(Tc)-1, length(Tc));
    for i = 1:length(Tc)-1
        Z(i, i:i+1) = [Tc(i) -Tc(i+1)]*sm;
    end

    M = [M; Z];
    Fr = [Fr, zeros(1, size(Z, 1))];
    
    wR = ones(1, numel(inclR(inclR)))*param.wFreq*param.wR;
    wL = ones(1, numel(inclL(inclL)))*param.wFreq*param.wL;
    wT = ones(1, size(Z, 1))*param.wTemp;
    w = [wR wL wT];
    
    x = lscov(M, Fr', w');
    %x(x <= 0.2 | x > 5) = 1;
    x = min(param.expMax, max(x, param.expMin));
    x = x';

    relax = param.relax;
    
%     rsum = resRrel + resLrel;
%     if rsum < 0.04
%         relax = (rsum * 25);
%     else
%         relax = relax0;
%     end
end

end

%--------------------------------------------------------------------------
function [h0, hc] = l_updatePlot(hAx, h0, hc, x0, y0, xc, yc, color)

if ~isempty(h0)
    delete(h0);
end
if ~isempty(hc)
    set(hc, 'color', [0.9 0.9 0.9])
end
hc = plot(hAx, xc, yc, 'Color', color, 'LineWidth', 2);
if ~isempty(h0)
    h0 = plot(hAx, x0, y0, 'Color', 'g');
end

end
