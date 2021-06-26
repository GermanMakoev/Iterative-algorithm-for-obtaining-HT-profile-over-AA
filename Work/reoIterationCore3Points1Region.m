function reoIterationCore3Points1Region(hLib, mfoData, Mmap, freqs, Robs, Lobs, H0, T0, Hb, Tb, Ht1, Ht2, Hc, Tc, NTv, posR, posL, param, Points)

nNT = 1;
NT = NTv(nNT);

outCntres = zeros(length(NTv));
outCntNT = zeros(length(NTv));
outRes = zeros(length(NTv));

[~, hAxes] = ssuCreateMultiAxesFig(3, 3, '', false, false);
hT0 = [];
title(hAxes(1,1),'Температурно-высотный профиль полутени');
xlabel(hAxes(1,1),'Высота, см');
ylabel(hAxes(1,1),'Температура, К');
title(hAxes(1,2),'Зависимость невязки от числа итераций');
xlabel(hAxes(1,2),'Число итераций');
ylabel(hAxes(1,2),'Невязка');
title(hAxes(2,1),'Спектры, правая поляризация, 1-я точка');
xlabel(hAxes(2,1),'Частота, Гц');
ylabel(hAxes(2,1),'Поток, с.е.п./угл.сек.');
title(hAxes(2,2),'Спектры, правая поляризация, 2-я точка');
xlabel(hAxes(2,2),'Частота, Гц');
ylabel(hAxes(2,2),'Поток, с.е.п./угл.сек.');
title(hAxes(2,3),'Спектры, правая поляризация, 3-я точка');
xlabel(hAxes(2,3),'Частота, Гц');
ylabel(hAxes(2,3),'Поток, с.е.п./угл.сек.');
title(hAxes(3,1),'Спектры, левая поляризация, 1-я точка');
xlabel(hAxes(3,1),'Частота, Гц');
ylabel(hAxes(3,1),'Поток, с.е.п./угл.сек.');
title(hAxes(3,2),'Спектры, левая поляризация, 2-я точка');
xlabel(hAxes(3,2),'Частота, Гц');
ylabel(hAxes(3,2),'Поток, с.е.п./угл.сек.');
title(hAxes(3,3),'Спектры, левая поляризация, 3-я точка');
xlabel(hAxes(3,3),'Частота, Гц');
ylabel(hAxes(3,3),'Поток, с.е.п./угл.сек.');
if ~isempty(H0) && ~isempty(T0)
    hT0 = plot(hAxes(1,1), H0, T0, 'Color', 'g');
end
set(hAxes(1,1), 'YScale', 'log', 'XLim', [0 1e9])
set(hAxes(1,2), 'YScale', 'log')
set(hAxes,'Fontsize',16)
hTc = [];
hRs0_1 = plot(hAxes(2,1), freqs, Robs(:,1), 'Color', 'g');
hRsc_1 = [];
hLs0_1 = plot(hAxes(3,1), freqs, Lobs(:,1), 'Color', 'g');
hLsc_1 = [];
hRs0_2 = plot(hAxes(2,2), freqs, Robs(:,2), 'Color', 'g');
hRsc_2 = [];
hLs0_2 = plot(hAxes(3,2), freqs, Lobs(:,2), 'Color', 'g');
hLsc_2 = [];
hRs0_3 = plot(hAxes(2,3), freqs, Robs(:,3), 'Color', 'g');
hRsc_3 = [];
hLs0_3 = plot(hAxes(3,3), freqs, Lobs(:,3), 'Color', 'g');
hLsc_3 = [];
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
    cntR = zeros(length(freqs),length(Hc), 3);
    cntL = zeros(length(freqs),length(Hc), 3);
    
    Tc = relax*x.*Tc + (1-relax)*Tprev;
    Tc(Tc < param.Tmin) = param.Tmin;
    Tprev = Tc;
    
    Hcalc = [Hb Hc];
    Tcalc = [Tb Tc];
    Dcalc = NT./Tcalc;
    
    R = zeros(length(freqs), length(Hc), 3);
    L = zeros(length(freqs), length(Hc), 3);
    HRm = zeros(length(freqs), length(Hc), 3);
    HLm = zeros(length(freqs), length(Hc), 3);
    
    Rcalc = zeros(length(freqs),3);
    Lcalc = zeros(length(freqs),3);
    posRcalc = zeros(length(freqs),3);
    posLcalc = zeros(length(freqs),3);
for l=1:3
    for k = 1:length(freqs)
            posRcalc(k,l) = posR(k,l);
            posLcalc(k,l) = posL(k,l);
        
        [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
         depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft, ...
         pScanRight, pScanLeft] = ...
                 reoCalculate(hLib, mfoData, Hcalc, Tcalc, Dcalc, Mmap, freqs(k), param.harms, param.pTauL, param.mode, param.c, param.b, posRcalc(k,l), posLcalc(k,l));
        Rcalc(k,l) = pScanRight(Points(l));
        Lcalc(k,l) = pScanLeft(Points(l));
                 
        for i= 1:size(pFRight, 1)
            for j = 1:size(pFRight, 2)
                for m = 1:length(Hc)
                    if depthRight(i,j) >= 1 && Ht1(m) <= pHLRight(i,j,1) && pHLRight(i,j,1) < Ht2(m)
                        R(k,m,l) = R(k,m,l) + pFLRight(i,j,1);
                        cntR(k,m,l) = cntR(k,m,l) + 1;
                        HRm(k,m,l) = HRm(k,m,l) + pHLRight(i,j,1);
                    end
                    if depthLeft(i,j) >= 1 && Ht1(m) <= pHLLeft(i,j,1) && pHLLeft(i,j,1) < Ht2(m)
                        L(k,m,l) = L(k,m,l) + pFLLeft(i,j,1);
                        cntL(k,m,l) = cntL(k,m,l) + 1;
                        HLm(k,m,l) = HLm(k,m,l) + pHLLeft(i,j,1);
                    end
                end
            end
        end
        HRm(k,:,l) = HRm(k,:,l)./cntR(k,:,l);
        HLm(k,:,l) = HLm(k,:,l)./cntL(k,:,l);
    end
end  

    devR=(Robs-Rcalc).^2;
    devL=(Lobs-Lcalc).^2;
    resR = sqrt(sum(sum(devR)));
    resRrel = resR/sum(sum(Robs));
    resL = sqrt(sum(sum(devL)));
    resLrel = resL/sum(sum(Lobs));
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
    [hRs0_1, hRsc_1] = l_updatePlot(hAxes(2, 1), hRs0_1, hRsc_1, freqs, Robs(:,1), freqs, Rcalc(:,1), 'r');
    [hLs0_1, hLsc_1] = l_updatePlot(hAxes(3, 1), hLs0_1, hLsc_1, freqs, Lobs(:,1), freqs, Lcalc(:,1), 'b');
    [hRs0_2, hRsc_2] = l_updatePlot(hAxes(2, 2), hRs0_2, hRsc_2, freqs, Robs(:,2), freqs, Rcalc(:,2), 'r');
    [hLs0_2, hLsc_2] = l_updatePlot(hAxes(3, 2), hLs0_2, hLsc_2, freqs, Lobs(:,2), freqs, Lcalc(:,2), 'b');
    [hRs0_3, hRsc_3] = l_updatePlot(hAxes(2, 3), hRs0_3, hRsc_3, freqs, Robs(:,3), freqs, Rcalc(:,3), 'r');
    [hLs0_3, hLsc_3] = l_updatePlot(hAxes(3, 3), hLs0_3, hLsc_3, freqs, Lobs(:,3), freqs, Lcalc(:,3), 'b');
    lr = length(resRv);
    if lr > 1
        plot(hAxes(1,2), [lr-1 lr], resRv(end-1:end), 'Color', 'r')
        plot(hAxes(1,2), [lr-1 lr], resLv(end-1:end), 'Color', 'b')
        plot(hAxes(1,2), [lr-1 lr], resRv(end-1:end)+resLv(end-1:end), 'Color', 'k')
    end
    drawnow
    
    R = [R(:, :, 1); R(:, :,  2); R(:, :,  3)];
    L = [L(:, :, 1); L(:, :,  2); L(:, :,  3)];
    Robs = [Robs(:,1); Robs(:,2); Robs(:,3)];
    Lobs = [Lobs(:,1); Lobs(:,2); Lobs(:,3)];
    Rcalc = [Rcalc(:,1); Rcalc(:,2); Rcalc(:,3)];
    Lcalc = [Lcalc(:,1); Lcalc(:,2); Lcalc(:,3)];
    
    inclR = true(1, 3*length(freqs));
    for i = 1:3*length(freqs)
        if all(R(i, :) == 0)
            inclR(i) = false;
        end
    end
    inclL = true(1, 3*length(freqs));
    for i = 1:3*length(freqs)
        if all(L(i, :) == 0)
            inclL(i) = false;
        end
    end
    
    M = [R(inclR, :); L(inclL, :)];
    Fobs =  [Robs(inclR);  Lobs(inclL)];
    Fcalc = [Rcalc(inclR); Lcalc(inclL)];
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
    %param.expMax=5;
    x = min(param.expMax, max(x, param.expMin));
    x = x';

    relax = param.relax;
    
    Robs=[Robs(1:length(freqs)),Robs(length(freqs)+1:2*length(freqs)),Robs(2*length(freqs)+1:3*length(freqs))];
    Lobs=[Lobs(1:length(freqs)),Lobs(length(freqs)+1:2*length(freqs)),Lobs(2*length(freqs)+1:3*length(freqs))];
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
