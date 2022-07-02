Tcalc8 = load('Tcalc8.mat','Tcalc');
Tcalc9 = load('Tcalc9.mat','Tcalc');
Tcalc10 = load('Tcalc10.mat','Tcalc');
Tcalc11 = load('Tcalc11.mat','Tcalc');
Tcalc12 = load('Tcalc12.mat','Tcalc');
%load('Temp0.mat','Temp0');
%load('H0.mat','H0');
load('Hcalc.mat','Hcalc');
%plot(H0,Temp0);
%hold on;
plot(Hcalc,Tcalc8.Tcalc);
hold on;
plot(Hcalc,Tcalc9.Tcalc);
hold on;
plot(Hcalc,Tcalc10.Tcalc);
hold on;
plot(Hcalc,Tcalc11.Tcalc);
hold on;
plot(Hcalc,Tcalc12.Tcalc);
title('Влияние магнитного поля на профиль')
xlabel('Высота,м')
ylabel('Температура,К')
%%
x = Hcalc;
tiledlayout(2,1)

% Top plot
ax1 = nexttile;
plot(ax1,Hcalc,Tcalc8.Tcalc(1,:));
hold on;
plot(ax1,Hcalc,Tcalc9.Tcalc(1,:));
hold on;
plot(ax1,Hcalc,Tcalc10.Tcalc(1,:));
hold on;
plot(ax1,Hcalc,Tcalc11.Tcalc(1,:));
hold on;
plot(ax1,Hcalc,Tcalc12.Tcalc(1,:));
title(ax1,'Влияние магнитного поля на профиль полутени')
xlabel(ax1,'Высота,м')
ylabel(ax1,'Температура,К')

% Bottom plot
ax2 = nexttile;
plot(ax2,Hcalc,Tcalc8.Tcalc(2,:));
hold on;
plot(ax2,Hcalc,Tcalc9.Tcalc(2,:));
hold on;
plot(ax2,Hcalc,Tcalc10.Tcalc(2,:));
hold on;
plot(ax2,Hcalc,Tcalc11.Tcalc(2,:));
hold on;
plot(ax2,Hcalc,Tcalc12.Tcalc(2,:));
title(ax2,'Влияние магнитного поля на профиль тени')
xlabel(ax1,'Высота,м')
ylabel(ax1,'Температура,К')
%%
NT = [6e14, 2e15, 6e15, 2e16, 6e16];
NB = [7.07, 3.74, 1.04, 4.52, 7.48];
plot(NT,NB);
title('Зависимость невязки от NT')
xlabel('NT, K*см^{-3}')
ylabel('Невязка,%')
%%
B = [0.8, 0.9, 1.0, 1.1, 1.2];
NB = [7.28, 2.20, 1.04, 1.38, 1.89];
plot(B,NB);
title('Зависимость невязки от магнитного поля')
xlabel('B/B_{0}')
ylabel('Невязка,%')
%%
Tcalc = load('T_12419.mat','Tcalc');
load('Hcalc.mat','Hcalc');
plot(Hcalc,Tcalc.Tcalc(1,:));
hold on;
plot(Hcalc,Tcalc.Tcalc(2,:));
title('ТВП, 12419')
xlabel('Высота, см')
ylabel('Поток, с.е.п./угл.с.')
