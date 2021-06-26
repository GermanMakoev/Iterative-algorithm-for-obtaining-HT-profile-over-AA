function reoCalcMapScanConveyorTest

mfoData = iouLoadField('C:\Users\Alan Makoev\Desktop\Matlab\dip_140_05_16e8.mat');

if isempty(mfoData)
    return
end

QT = 1;
freefree = 0;

hLib = reoInitLibrary(QT, freefree);
if hLib == 0
    return
end

gridstep = [1 1]; % in arcsec
posangle = 0;
vcos = [0; 0; 1];

[M, base] = reoSetField(hLib, mfoData, gridstep, posangle, vcos);

%Исходная модель
format short e
H = [1, 1e8, 1.1e8, 2e10];
Temp = [1e4, 1e4, 5e6, 5e6];
Dens = 3e15./Temp;
freqs = [3:0.5:18]*1e9;
harms = 2:4;
pTauL = [1 25];
mode = 2;
c = 0;
b = 0;

R_real = zeros(1,length(freqs));
L_real = zeros(1,length(freqs));
R_pos = zeros(1,length(freqs));
L_pos = zeros(1,length(freqs));
for k = 1:length(freqs)
    [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
     depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft, ...
     pScanRight, pScanLeft] = ...
             reoCalculate(hLib, mfoData, H, Temp, Dens, M, freqs(k), harms, pTauL, mode, c, b);
    [R_real(k), R_pos(k)] = max(pScanRight);
    [L_real(k), L_pos(k)] = max(pScanLeft);
end

%Модель для варьирования
Ht = [1 1e8 1.25e8 1.5e8 1.75e8 2e8 2.25e8 2.5e8 3e8 6e8 10e8 15e8 20e8 2e10];
Hc = (Ht(1:end-1)+Ht(2:end))/2;
Tc = 1e4*ones(1, length(Hc));
Tc(3:end) = 3e6;
Dc = 3e15./Tc;
x=ones(1,length(Tc));

%Цикл итераций
M_1=2;
N=zeros(1,M_1);
for z=1:M_1
    for h=1:length(x)
    Tc(h)=x(h)*Tc(h);
    end
    R = zeros(length(freqs),length(Hc));
    L = zeros(length(freqs),length(Hc));
%Полные потоки  
    for k = 1:length(freqs)
        [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
         depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft, ...
         pScanRight, pScanLeft] = ...
                 reoCalculate(hLib, mfoData, Hc, Tc, Dc, M, freqs(k), harms, pTauL, mode, c, b);
        [R_full(k), posR] = max(pScanRight);
        [L_full(k), posL] = max(pScanLeft);
        [depthRight, pFRight, pTauRight, pHLRight, pFLRight, psLRight, ...
         depthLeft, pFLeft, pTauLeft, pHLLeft, pFLLeft, psLLeft, ...
         pScanRight, pScanLeft] = ...
                 reoCalculate(hLib, mfoData, Hc, Tc, Dc, M, freqs(k), harms, pTauL, mode, c, b, posR, posL);
%Вклады от оптически толстых слоёв        
      for i= 1:size(pFRight, 1)
            for j = 1:size(pFRight, 2)
                for m = 1:length(Hc)
                    if depthRight(i,j) >= 1 && Ht(m) <= pHLRight(i,j,1) && pHLRight(i,j,1) < Ht(m+1)
                        R(k,m) = R(k,m) + pFLRight(i,j,1);
                    end
                    if depthLeft(i,j) >= 1 && Ht(m) <= pHLLeft(i,j,1) && pHLLeft(i,j,1) < Ht(m+1)
                        L(k,m) = L(k,m) + pFLLeft(i,j,1);
                    end
                end
            end
        end
    end  
    R;
    L;
%Вырезаем нулевые строки    
    R(18:end,:) = [];
    L(14:end,:) = [];
    R_real(18:end) = [];
    L_real(14:end) = [];
    R_full(18:end) = [];
    L_full(14:end) = [];
    F=[R;L];
    size(F);
    F_real=[R_real L_real];
    F_full=[R_full L_full];
%Вычисляем правую часть
    for i = 1:length(F_real)
        F_1(i)=F_real(i)+sum(F(i,:))-F_full(i);
    end
%Добавляем уравнения для температур
    for i=1:length(Tc)
           Z(i,i)=Tc(i)*1e-9;
    end
    for i=1:(length(Tc)-1)
           Z(i,i+1)=-Tc(i+1)*1e-9;
    end
    Z;
    Z(13,:)=[];
    F_2=[F;Z]
%Добавляем веса
    w=ones(1,length(F_2));
    w(31:end)=0.2;
    %w(15:17)=0.1;
    %w(28:30)=0.1;
%Решаем систему
    F_3=[F_1 zeros(1,12)];
 x=lscov(F_2,F_3',w);
    for i=1:length(x)
        if x(i)<0
            x(i)=abs(x(i));
            %x(i)=1;
        end
    end
    x
%Вычисляем невязку
B=F_2*x;
for i=1:length(x)
    S(i)=B(i)-F_3(i);
end
S
N(z)=sum(w(i)*(S(i))^2);
N(z)
%Выводим температуры и номер итерации 
 Tc
 z
end
 subplot(2,1,1), plot(N)
 subplot(2,1,2), plot(freqs(1:17),R_full)
 hold on;
 plot(freqs(1:17),R_real)
 legend('R_{full}','R_{real}')
 utilsFreeLibrary(hLib);
end

