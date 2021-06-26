function reoIterationCore1Point2RegionV2(hLib, mfoData, Mmap, freqs, Robs, Lobs, H0, T0, Hb, Tb, Ht1, Ht2, Hc, Tc, NTv, posR, posL, Mask, param, diagrH, diagrV, gridstep)

nNT = 1;
NT = NTv(nNT);
nModels=2;

outCntres = zeros(length(NTv));
outCntNT = zeros(length(NTv));
outRes = zeros(length(NTv));

[~, hAxes] = ssuCreateMultiAxesFig(2, 3, '', false, false);
hT0_1 = [];
hT0_2 = [];
title(hAxes(1,1),'ТВП полутени');
xlabel(hAxes(1,1),'Высота, см');
ylabel(hAxes(1,1),'Температура, К');
title(hAxes(1,2),'ТВП тени');
xlabel(hAxes(1,2),'Высота, см');
ylabel(hAxes(1,2),'Температура, К');
title(hAxes(1,3),'Зависимость невязки от числа итераций');
xlabel(hAxes(1,3),'Число итераций');
ylabel(hAxes(1,3),'Невязка');
title(hAxes(2,1),'Спектры, правая поляризация');
xlabel(hAxes(2,1),'Частота, Гц');
ylabel(hAxes(2,1),'Поток, с.е.п./угл.сек.');
title(hAxes(2,2),'Спектры, левая поляризация');
xlabel(hAxes(2,2),'Частота, Гц');
ylabel(hAxes(2,2),'Поток, с.е.п./угл.сек.');
%title(hAxes(1,1),'Height-temperature profile');
%xlabel(hAxes(1,1),'Height, cm');
%ylabel(hAxes(1,1),'Temperature, K');
%title(hAxes(1,2),'Height-temperature profile');
%xlabel(hAxes(1,2),'Height, cm');
%ylabel(hAxes(1,2),'Temperature, K');
%title(hAxes(1,3),'Dependence of the deviation on the number of iterations');
%xlabel(hAxes(1,3),'Number of iterations');
%ylabel(hAxes(1,3),'Deviation');
%title(hAxes(2,1),'Spectra, right polarization');
%xlabel(hAxes(2,1),'Frequency, Hz');
%ylabel(hAxes(2,1),'Flux, s.f.u./arc.sec.');
%title(hAxes(2,2),'Spectra, left polarization');
%xlabel(hAxes(2,2),'Frequency, Hz');
%ylabel(hAxes(2,2),'Flux, s.f.u./arc.sec.');
if ~isempty(H0) && ~isempty(T0(1,:))
    hT0_1 = plot(hAxes(1,1), H0, T0(1,:), 'Color', 'g');
end
if ~isempty(H0) && ~isempty(T0(2,:))
    hT0_2 = plot(hAxes(1,2), H0, T0(2,:), 'Color', 'g');
end
set(hAxes(1,1), 'YScale', 'log', 'XLim', [0 1e9])
set(hAxes(1,2), 'YScale', 'log', 'XLim', [0 1e9])
set(hAxes(1,3), 'YScale', 'log')
set(hAxes,'Fontsize',16)
hTc_1 = [];
hTc_2 = [];
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
x = ones(2, length(Hc));
while true    
    cntR = zeros(length(freqs), length(Hc), nModels);
    cntL = zeros(length(freqs), length(Hc), nModels);
    
    Tc(1,:) = relax*x(1,:).*Tc(1,:) + (1-relax)*Tprev(1,:);
    Tc(2,:) = relax*x(2,:).*Tc(2,:) + (1-relax)*Tprev(2,:);
    Tc(Tc < param.Tmin) = param.Tmin;
    Tprev = Tc;
    
    Hcalc = [Hb Hc];
    Tcalc(1,:) = [Tb(1,:) Tc(1,:)];
    Tcalc(2,:) = [Tb(2,:) Tc(2,:)];
    Dcalc = NT./Tcalc;
    
    R = zeros(length(freqs), length(Hc), nModels);
    L = zeros(length(freqs), length(Hc), nModels);
    HRm = zeros(length(freqs), length(Hc), nModels);
    HLm = zeros(length(freqs), length(Hc), nModels);
    
    Rcalc = zeros(1, length(freqs));
    Lcalc = zeros(1, length(freqs));
    for k = 1:length(freqs)
        depthRightW=zeros(Mmap);
        pFRightW = zeros(Mmap);
        pFLRightW = zeros(Mmap);
        pHLRightW=zeros(Mmap);
        depthLeftW=zeros(Mmap);
        pFLeftW = zeros(Mmap);
        pFLLeftW = zeros(Mmap);
        pHLLeftW=zeros(Mmap);
%        if isempty(posL) || isempty(posR)
        for m = 1:nModels
            [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
            depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft] = ...
                     reoCalculate(hLib, mfoData, Hcalc, Tcalc(m,:), Dcalc(m,:), Mmap, freqs(k), param.harms, param.pTauL, param.mode, param.c, param.b);
            depthRightW(Mask==m) = depthRight(Mask==m); 
            pFRightW(Mask == m) = pFRight(Mask == m);
            pR = pFLRight(:, :, 1);
            pFLRightW(Mask == m) = pR(Mask == m);
            pHR=pHLRight(:, :, 1);
            pHLRightW(Mask == m) = pHR(Mask == m);
            depthLeftW(Mask==m) = depthLeft(Mask==m); 
            pFLeftW(Mask == m) = pFLeft(Mask == m);
            pL = pFLLeft(:, :, 1);
            pFLLeftW(Mask == m) = pL(Mask == m);
            pHL=pHLLeft(:, :, 1);
            pHLLeftW(Mask == m) = pHL(Mask == m); 
        end     
%             scanR = gstMapConv(pFRightW, diagrH(k,:), diagrV(k,:), gridstep);
%             scanL = gstMapConv(pFLeftW, diagrH(k,:), diagrV(k,:), gridstep);
%             [~, posRcalc] = max(scanR);
%             
%             [~, posLcalc] = max(scanL);
%         else
%             posRcalc = posR(k);
%             posLcalc = posL(k);  
%         end    
         
        posRcalc = posR(k);
        posLcalc = posL(k);  
        pFLRightW = gstMapMult(pFLRightW, diagrH(k,:), diagrV(k,:), gridstep, posRcalc);
        pFLLeftW = gstMapMult(pFLLeftW, diagrH(k,:), diagrV(k,:), gridstep, posLcalc);
        scanR = gstMapConv(pFRightW, diagrH(k,:), diagrV(k,:), gridstep);
        scanL = gstMapConv(pFLeftW, diagrH(k,:), diagrV(k,:), gridstep);
        Rcalc(k) = max(scanR);
        Lcalc(k) = max(scanL);
                 
        for i= 1:size(pFRightW, 1)
            for j = 1:size(pFRightW, 2)
                if Mask(i,j)>0
                u = Mask(i, j);
                for m = 1:length(Hc)
                    if depthRightW(i,j) >= 1 && Ht1(m) <= pHLRightW(i,j) && pHLRightW(i,j) < Ht2(m)
                        R(k,m,u) = R(k,m,u) + pFLRightW(i,j);
                        cntR(k,m,u) = cntR(k,m,u) + 1;
                        HRm(k,m,u) = HRm(k,m,u)+ pHLRightW(i,j);
                    end
                    if depthLeftW(i,j) >= 1 && Ht1(m) <= pHLLeftW(i,j) && pHLLeftW(i,j) < Ht2(m)
                        L(k,m,u) = L(k,m,u) + pFLLeftW(i,j);
                        cntL(k,m,u) = cntL(k,m,u) + 1;
                        HLm(k,m,u) = HLm(k,m,u) + pHLLeftW(i,j);
                    end
                end
                end
            end
        end
        HRm(k,:) = HRm(k,:)./cntR(k,:);
        HLm(k,:) = HLm(k,:)./cntL(k,:);
    end
    
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
    
    if isempty(T0)
        [hT0_1, hTc_1] = l_updatePlot(hAxes(1, 1), hT0_1, hTc_1, H0, [], Hcalc, Tcalc(1,:), 'm');
        [hT0_2, hTc_2] = l_updatePlot(hAxes(1, 2), hT0_2, hTc_2, H0, [], Hcalc, Tcalc(2,:), 'm');
    else
        [hT0_1, hTc_1] = l_updatePlot(hAxes(1, 1), hT0_1, hTc_1, H0, T0(1,:), Hcalc, Tcalc(1,:), 'm');
        [hT0_2, hTc_2] = l_updatePlot(hAxes(1, 2), hT0_2, hTc_2, H0, T0(2,:), Hcalc, Tcalc(2,:), 'm');
    end
    [hRs0, hRsc] = l_updatePlot(hAxes(2, 1), hRs0, hRsc, freqs, Robs, freqs, Rcalc, 'r');
    [hLs0, hLsc] = l_updatePlot(hAxes(2, 2), hLs0, hLsc, freqs, Lobs, freqs, Lcalc, 'b');
    lr = length(resRv);
    if lr > 1
        plot(hAxes(1,3), [lr-1 lr], resRv(end-1:end), 'Color', 'r')
        plot(hAxes(1,3), [lr-1 lr], resLv(end-1:end), 'Color', 'b')
        plot(hAxes(1,3), [lr-1 lr], resRv(end-1:end)+resLv(end-1:end), 'Color', 'k')
    end
    drawnow
    
    R = [R(:, :, 1) R(:, :, 2)];
    L = [L(:, :, 1) L(:, :, 2)];
    
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
    
    Z1 = zeros(length(Tc(1,:))-1, length(Tc(1,:)));
    for i = 1:length(Tc(1,:))-1
        Z1(i,i:i+1) = [Tc(1,i) -Tc(1,i+1)]*sm;
    end
 
    Z2 = zeros(length(Tc(2,:))-1, length(Tc(2,:)));
    for i = 1:length(Tc(2,:))-1
        Z2(i, i:i+1) = [Tc(2,i) -Tc(2,i+1)]*sm;
    end
    
    D1=zeros(size(Tc(1,:)));
    for i=1:length(Tc(1,:))
    D1(i,i)=Tc(1,i)*sm;
    end
    
    D2=zeros(size(Tc(2,:)));
    for i=1:length(Tc(2,:))
    D2(i,i)=-Tc(2,i)*sm;
    end
    
    D=[D1,D2];
    
    Z=[Z1,zeros(size(Z2));zeros(size(Z1)), Z2];

    M = [M ; Z; D];
    Fr = [Fr, zeros(1, size(Z, 1)),zeros(1, size(D, 1))];
    
    wR = ones(1, numel(inclR(inclR)))*param.wFreq*param.wR;
    wL = ones(1, numel(inclL(inclL)))*param.wFreq*param.wL;
    wT = ones(1, size(Z, 1))*param.wTemp;
    wT_D = ones(1, size(D, 1))*param.wTemp_D;
    w = [wR wL wT wT_D];
    
    t = lscov(M, Fr', w');
    t = min(param.expMax, max(t, param.expMin));
    x(1, :) = t(1:length(Hc));
    x(2, :) = t(length(Hc) + (1:length(Hc)));
%     x(1,:) = lscov(M(:,1:length(Hc)), Fr', w');
%     x(2,:) = lscov(M(:,length(Hc)+1:2*length(Hc)), Fr', w');
%     %x(x <= 0.2 | x > 5) = 1;
%     x(1,:) = min(param.expMax, max(x(1,:), param.expMin));
%     x(2,:) = min(param.expMax, max(x(2,:), param.expMin));
%     x(1,:) = x(1,:)';
%     x(2,:) = x(2,:)';

    relax = param.relax;
end
    
%     rsum = resRrel + resLrel;
%     if rsum < 0.04
%         relax = (rsum * 25);
%     else
%         relax = relax0;
%     end
save('Tcalc8.mat','Tcalc');
save('Hcalc.mat','Hcalc');
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

