function reoIterationCore3Points2RegionsV2(hLib, mfoData, Mmap, freqs, freqsR, freqsL, Robs, Lobs, RB, LB, H0, T0, Hb, Tb, Ht1, Ht2, Hc, Tc, NTv, posR, posL, Mask, param, diagrH, diagrV, gridstep, Points)

nNT = 1;
NT = NTv(nNT);
nModels=2;

outCntres = zeros(length(NTv));
outCntNT = zeros(length(NTv));
outRes = zeros(length(NTv));

[~, hAxes] = ssuCreateMultiAxesFig(3, 3, '', false, false);
hT0_1 = [];
hT0_2 = [];
title(hAxes(1,1),'??? ????????');
xlabel(hAxes(1,1),'??????, ??');
ylabel(hAxes(1,1),'???????????, ?');
title(hAxes(1,2),'??? ????');
xlabel(hAxes(1,2),'??????, ??');
ylabel(hAxes(1,2),'???????????, ?');
title(hAxes(1,3),'??????? (????? ????????)');
xlabel(hAxes(1,3),'????? ????????');
ylabel(hAxes(1,3),'???????');
title(hAxes(2,1),'???????, ??, 1-? ?????');
xlabel(hAxes(2,1),'???????, ??');
ylabel(hAxes(2,1),'?????, ?.?.?./???.???.');
title(hAxes(2,2),'???????, ??, 2-? ?????');
xlabel(hAxes(2,2),'???????, ??');
ylabel(hAxes(2,2),'?????, ?.?.?./???.???.');
title(hAxes(2,3),'???????, ??, 3-? ?????');
xlabel(hAxes(2,3),'???????, ??');
ylabel(hAxes(2,3),'?????, ?.?.?./???.???.');
title(hAxes(3,1),'???????, ??, 1-? ?????');
xlabel(hAxes(3,1),'???????, ??');
ylabel(hAxes(3,1),'?????, ?.?.?./???.???.');
title(hAxes(3,2),'???????, ??, 2-? ?????');
xlabel(hAxes(3,2),'???????, ??');
ylabel(hAxes(3,2),'?????, ?.?.?./???.???.');
title(hAxes(3,3),'???????, ??, 3-? ?????');
xlabel(hAxes(3,3),'???????, ??');
ylabel(hAxes(3,3),'?????, ?.?.?./???.???.');
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
hRs0_1 = plot(hAxes(2,1), freqsR, Robs(:,1), 'Color', 'g');
hRsc_1 = [];
hLs0_1 = plot(hAxes(3,1), freqsL, Lobs(:,1), 'Color', 'g');
hLsc_1 = [];
hRs0_2 = plot(hAxes(2,2), freqsR, Robs(:,2), 'Color', 'g');
hRsc_2 = [];
hLs0_2 = plot(hAxes(3,2), freqsL, Lobs(:,2), 'Color', 'g');
hLsc_2 = [];
hRs0_3 = plot(hAxes(2,3), freqsR, Robs(:,3), 'Color', 'g');
hRsc_3 = [];
hLs0_3 = plot(hAxes(3,3), freqsL, Lobs(:,3), 'Color', 'g');
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
    Rcalc = zeros(length(freqs), 3);
    Lcalc = zeros(length(freqs), 3);
    posRcalc = zeros(length(freqs), 3);
    posLcalc = zeros(length(freqs), 3);
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
    scanRs = [];
    scanLs = [];
    for k = 1:length(freqs)
        depthRightW = zeros(Mmap);
        pFRightW = zeros(Mmap);
        pFLRightW = zeros(Mmap);
        pHLRightW = zeros(Mmap);
        depthLeftW = zeros(Mmap);
        pFLeftW = zeros(Mmap);
        pFLLeftW = zeros(Mmap);
        pHLLeftW = zeros(Mmap);
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
         
        posRcalc(k,l) = posR(k,l);
        posLcalc(k,l) = posL(k,l);  
        pFLRightW = gstMapMult(pFLRightW, diagrH(k,:), diagrV(k,:), gridstep, posRcalc(k,l));
        pFLLeftW = gstMapMult(pFLLeftW, diagrH(k,:), diagrV(k,:), gridstep, posLcalc(k,l));
        scanR = gstMapConv(pFRightW, diagrH(k,:), diagrV(k,:), gridstep);
        scanL = gstMapConv(pFLeftW, diagrH(k,:), diagrV(k,:), gridstep);
        scanRs = [scanRs;scanR];
        scanLs = [scanLs;scanL];
        Rcalc(k,l) = scanR(Points(1,l));
        Lcalc(k,l) = scanL(Points(2,l));
                 
        for i= 1:size(pFRightW, 1)
            for j = 1:size(pFRightW, 2)
                if Mask(i,j)>0
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
        end
        HRm(k,:,l) = HRm(k,:,l)./cntR(k,:,l);
        HLm(k,:,l) = HLm(k,:,l)./cntL(k,:,l);
    end
    save(['scanR.mat'],'scanRs');
    save(['scanL.mat'],'scanLs');
end

Rcalc((length(freqsR)+1):end,:) = [];
R((length(freqsR)+1):end,:,:,:) = [];
Lcalc((length(freqsL)+1):end,:) = [];
L((length(freqsL)+1):end,:,:,:) = []; 

Rcalc = Rcalc + RB;
Lcalc = Lcalc + LB;
                    
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
    
    if isempty(T0)
        [hT0_1, hTc_1] = l_updatePlot(hAxes(1, 1), hT0_1, hTc_1, H0, [], Hcalc, Tcalc(1,:), 'm');
        [hT0_2, hTc_2] = l_updatePlot(hAxes(1, 2), hT0_2, hTc_2, H0, [], Hcalc, Tcalc(2,:), 'm');
    else
        [hT0_1, hTc_1] = l_updatePlot(hAxes(1, 1), hT0_1, hTc_1, H0, T0(1,:), Hcalc, Tcalc(1,:), 'm');
        [hT0_2, hTc_2] = l_updatePlot(hAxes(1, 2), hT0_2, hTc_2, H0, T0(2,:), Hcalc, Tcalc(2,:), 'm');
    end
    [hRs0_1, hRsc_1] = l_updatePlot(hAxes(2, 1), hRs0_1, hRsc_1, freqsR, Robs(:,1), freqsR, Rcalc(:,1), 'r');
    [hLs0_1, hLsc_1] = l_updatePlot(hAxes(3, 1), hLs0_1, hLsc_1, freqsL, Lobs(:,1), freqsL, Lcalc(:,1), 'b');
    [hRs0_2, hRsc_2] = l_updatePlot(hAxes(2, 2), hRs0_2, hRsc_2, freqsR, Robs(:,2), freqsR, Rcalc(:,2), 'r');
    [hLs0_2, hLsc_2] = l_updatePlot(hAxes(3, 2), hLs0_2, hLsc_2, freqsL, Lobs(:,2), freqsL, Lcalc(:,2), 'b');
    [hRs0_3, hRsc_3] = l_updatePlot(hAxes(2, 3), hRs0_3, hRsc_3, freqsR, Robs(:,3), freqsR, Rcalc(:,3), 'r');
    [hLs0_3, hLsc_3] = l_updatePlot(hAxes(3, 3), hLs0_3, hLsc_3, freqsL, Lobs(:,3), freqsL, Lcalc(:,3), 'b');
    lr = length(resRv);
    if lr > 1
        plot(hAxes(1,3), [lr-1 lr], resRv(end-1:end), 'Color', 'r')
        plot(hAxes(1,3), [lr-1 lr], resLv(end-1:end), 'Color', 'b')
        plot(hAxes(1,3), [lr-1 lr], resRv(end-1:end)+resLv(end-1:end), 'Color', 'k')
    end
    drawnow
    
    R = [R(:, :, 1, 1) R(:, :, 2, 1); R(:, :, 1, 2) R(:, :, 2, 2); R(:, :, 1, 3) R(:, :, 2, 3)];
    L = [L(:, :, 1, 1) L(:, :, 2, 1); L(:, :, 1, 2) L(:, :, 2, 2); L(:, :, 1, 3) L(:, :, 2, 3)];
    Robs = [Robs(:,1); Robs(:,2); Robs(:,3)];
    Lobs = [Lobs(:,1); Lobs(:,2); Lobs(:,3)];
    Rcalc = [Rcalc(:,1); Rcalc(:,2); Rcalc(:,3)];
    Lcalc = [Lcalc(:,1); Lcalc(:,2); Lcalc(:,3)];
    
    inclR = true(1, 3*length(freqsR));
    for i = 1:3*length(freqsR)
        if all(R(i, :) == 0)
            inclR(i) = false;
        end
    end
    inclL = true(1, 3*length(freqsL));
    for i = 1:3*length(freqsL)
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
    
    Robs=[Robs(1:length(freqsR)),Robs(length(freqsR)+1:2*length(freqsR)),Robs(2*length(freqsR)+1:3*length(freqsR))];
    Lobs=[Lobs(1:length(freqsL)),Lobs(length(freqsL)+1:2*length(freqsL)),Lobs(2*length(freqsL)+1:3*length(freqsL))];
    %Robs=[Robs(1:71),Robs(72:142),Robs(143:213)];
    %Lobs=[Lobs(1:71),Lobs(72:142),Lobs(143:213)];
end
    
%     rsum = resRrel + resLrel;
%     if rsum < 0.04
%         relax = (rsum * 25);
%     else
%         relax = relax0;
%     end
%save('T_12419.mat','Tcalc');
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

