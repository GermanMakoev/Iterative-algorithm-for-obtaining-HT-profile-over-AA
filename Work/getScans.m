function getScans(hLib, mfoData, Mmap, freqs, freqsR, freqsL, Robs, Lobs, RB, LB, H0, T0, Hb, Tb, Ht1, Ht2, Hc, Tc, NTv, posR, posL, Mask, param, diagrH, diagrV, gridstep)

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
end




