function param = reoGetParam

param.wFreq = 1; 
param.wTemp = 10; 
param.wTemp_D = 10;
param.wR = 1; 
param.wL = 1; 
param.relax = 1; 

param.harms = 2:4;
param.mode = 3;
param.c = 0;
param.b = 0;
param.pTauL = [1 25];
param.Tmin = 6000;

param.reslim = 5e-4;
param.rescntstab = 10;
param.rescntmax = 30;
param.expMax = 1.5;
param.expMin = 1./param.expMax;

end
