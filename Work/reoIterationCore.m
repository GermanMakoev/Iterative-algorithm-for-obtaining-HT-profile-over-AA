function reoIterationCore(hLib, mfoData, Mmap, freqs, Robs, Lobs, H0, T0, Hb, Tb, Ht1, Ht2, Hc, Tc, NTv, posR, posL, param, diagrH, diagrV, gridstep)

nNT = 1;
NT = NTv(nNT);
[M, base, Bc] = reoSetField(hLib, mfoData, gridstep);
Mask = reoGetMaskA(sqrt(Bc.x.^2 + Bc.y.^2 + Bc.z.^2));
nModels=2;

outCntres = zeros(length(NTv));
outCntNT = zeros(length(NTv));
outRes = zeros(length(NTv));

[~, hAxes] = ssuCreateMultiAxesFig(3, 3, '', false, false);
hT0_1 = [];
hT0_2 = [];
title(hAxes(1,1),'Температурно-высотный профиль полутени');
xlabel(hAxes(1,1),'Высота, см');
ylabel(hAxes(1,1),'Температура, К');
title(hAxes(1,2),'Температурно-высотный профиль тени');
xlabel(hAxes(1,2),'Высота, см');
ylabel(hAxes(1,2),'Температура, К');
title(hAxes(1,3),'Зависимость невязки от числа итераций');
xlabel(hAxes(1,3),'Число итераций');
ylabel(hAxes(1,3),'Невязка');
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
hRs0_1 = plot(hAxes(2,1), freqs, Robs(1,:), 'Color', 'g');
hRsc_1 = [];
hLs0_1 = plot(hAxes(3,1), freqs, Lobs(1,:), 'Color', 'g');
hLsc_1 = [];
hRs0_2 = plot(hAxes(2,2), freqs, Robs(2,:), 'Color', 'g');
hRsc_2 = [];
hLs0_2 = plot(hAxes(3,2), freqs, Lobs(2,:), 'Color', 'g');
hLsc_2 = [];
hRs0_3 = plot(hAxes(2,3), freqs, Robs(3,:), 'Color', 'g');
hRsc_3 = [];
hLs0_3 = plot(hAxes(3,3), freqs, Lobs(3,:), 'Color', 'g');
hLsc_3 = [];
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
    Rcalc = zeros(3, length(freqs));
    Lcalc = zeros(3, length(freqs));
    cntR = zeros(length(freqs), length(Hc), nModels, 3);
    cntL = zeros(length(freqs), length(Hc), nModels, 3);
    
    Tc(1,:) = relax*x(1,:).*Tc(1,:) + (1-relax)*Tprev(1,:);
    Tc(2,:) = relax*x(2,:).*Tc(2,:) + (1-relax)*Tprev(2,:);
    Tc(Tc < param.Tmin) = param.Tmin;
    Tprev = Tc;
    
    Hcalc = [Hb Hc];
    Tcalc(1,:) = [Tb(1,:) Tc(1,:)];
    Tcalc(2,:) = [Tb(2,:) Tc(2,:)];
    Dcalc = NT./Tcalc;
    
    R = zeros(length(freqs), length(Hc), nModels, 3);
    L = zeros(length(freqs), length(Hc), nModels, 3);
    HRm = zeros(length(freqs), length(Hc), nModels, 3);
    HLm = zeros(length(freqs), length(Hc), nModels, 3);
    
for l=1:3
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
         
        posRcalc(l) = posR(l,k);
        posLcalc(l) = posL(l,k);  
        pFLRightW = gstMapMult(pFLRightW, diagrH(k,:), diagrV(k,:), gridstep, posRcalc(l));
        pFLLeftW = gstMapMult(pFLLeftW, diagrH(k,:), diagrV(k,:), gridstep, posLcalc(l));
        scanR = gstMapConv(pFRightW, diagrH(k,:), diagrV(k,:), gridstep);
        scanL = gstMapConv(pFLeftW, diagrH(k,:), diagrV(k,:), gridstep);
        %Rcalc(k) = max(scanR);
        %Lcalc(k) = max(scanL);
        Rcalc(l,k) = scanR(24+14*(l-2));
        Lcalc(l,k) = scanL(24+14*(l-2));
                 
        for i= 1:size(pFRightW, 1)
            for j = 1:size(pFRightW, 2)
                u = Mask(i, j);
                for m = 1:length(Hc)
                    if depthRightW(i,j) >= 1 && Ht1(m) <= pHLRightW(i,j) && pHLRightW(i,j) < Ht2(m)
                        R(k,m,u,l) = R(k,m,u,l) + pFLRightW(i,j);
                        cntR(k,m,u,l) = cntR(k,m,u,l) + 1;
                        HRm(k,m,u,l) = HRm(k,m,u,l)+ pHLRightW(i,j);
                    end
                    if depthLeftW(i,j) >= 1 && Ht1(m) <= pHLLeftW(i,j) && pHLLeftW(i,j) < Ht2(m)
                        L(k,m,u,l) = L(k,m,u,l) + pFLLeftW(i,j);
                        cntL(k,m,u,l) = cntL(k,m,u,l) + 1;
                        HLm(k,m,u,l) = HLm(k,m,u,l) + pHLLeftW(i,j);
                    end
                end
            end
        end
        HRm(k,:,l) = HRm(k,:,l)./cntR(k,:,l);
        HLm(k,:,l) = HLm(k,:,l)./cntL(k,:,l);
    end
end
   
    resR = sqrt(sum(sum((Robs-Rcalc).^2)));
    resRrel = resR/sum(sum(Robs));
    resL = sqrt(sum(sum((Lobs-Lcalc).^2)));
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
    
    if isempty(T0)
        [hT0_1, hTc_1] = l_updatePlot(hAxes(1, 1), hT0_1, hTc_1, H0, [], Hcalc, Tcalc(1,:), 'm');
        [hT0_2, hTc_2] = l_updatePlot(hAxes(1, 2), hT0_2, hTc_2, H0, [], Hcalc, Tcalc(2,:), 'm');
    else
        [hT0_1, hTc_1] = l_updatePlot(hAxes(1, 1), hT0_1, hTc_1, H0, T0(1,:), Hcalc, Tcalc(1,:), 'm');
        [hT0_2, hTc_2] = l_updatePlot(hAxes(1, 2), hT0_2, hTc_2, H0, T0(2,:), Hcalc, Tcalc(2,:), 'm');
    end
    [hRs0_1, hRsc_1] = l_updatePlot(hAxes(2, 1), hRs0_1, hRsc_1, freqs, Robs(1,:), freqs, Rcalc(1,:), 'r');
    [hLs0_1, hLsc_1] = l_updatePlot(hAxes(3, 1), hLs0_1, hLsc_1, freqs, Lobs(1,:), freqs, Lcalc(1,:), 'b');
    [hRs0_2, hRsc_2] = l_updatePlot(hAxes(2, 2), hRs0_2, hRsc_2, freqs, Robs(2,:), freqs, Rcalc(2,:), 'r');
    [hLs0_2, hLsc_2] = l_updatePlot(hAxes(3, 2), hLs0_2, hLsc_2, freqs, Lobs(2,:), freqs, Lcalc(2,:), 'b');
    [hRs0_3, hRsc_3] = l_updatePlot(hAxes(2, 3), hRs0_3, hRsc_3, freqs, Robs(3,:), freqs, Rcalc(3,:), 'r');
    [hLs0_3, hLsc_3] = l_updatePlot(hAxes(3, 3), hLs0_3, hLsc_3, freqs, Lobs(3,:), freqs, Lcalc(3,:), 'b');
    lr = length(resRv);
    if lr > 1
        plot(hAxes(1,3), [lr-1 lr], resRv(end-1:end), 'Color', 'r')
        plot(hAxes(1,3), [lr-1 lr], resLv(end-1:end), 'Color', 'b')
        plot(hAxes(1,3), [lr-1 lr], resRv(end-1:end)+resLv(end-1:end), 'Color', 'k')
    end
    drawnow
    
    R = [R(:, :, 1, 1) R(:, :, 2, 1); R(:, :, 1, 2) R(:, :, 2, 2); R(:, :, 1, 3) R(:, :, 2, 3);];
    L = [L(:, :, 1, 1) L(:, :, 2, 1); L(:, :, 1, 2) L(:, :, 2, 2); L(:, :, 1, 3) L(:, :, 2, 3);];
    
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
    D1(i,i)=Tc(1,i);
    end
    
    D2=zeros(size(Tc(2,:)));
    for i=1:length(Tc(2,:))
    D2(i,i)=-Tc(2,i);
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

