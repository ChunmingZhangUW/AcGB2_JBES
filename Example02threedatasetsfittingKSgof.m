
clear all;

rand('seed',5252921);
burnin = 100; nt = 10000;
Qsp100 = xlsread('SP100_processed16.csv', 'SP100_processed16', 'C2:C1258');
nsp100 = length(Qsp100); Qsp100 = exp(Qsp100);

thetasp100acf = xlsread('ThreeDatasetsparameterestimations.xlsx','AcF','F1:F9');
[alpha,sigma,Frfit] = fn_recover_AcF(thetasp100acf(1:4),thetasp100acf(5:8),thetasp100acf(9),Qsp100,nsp100);
for i=1:nt
Fr = 1./exprnd(1,nsp100,1);
[h,pFrsp100(i)] = kstest2(Fr(burnin:nsp100),Frfit(burnin:nsp100),'Alpha',0.01);
end

thetasp100gb2 = xlsread('ThreeDatasetsparameterestimations.xlsx','AcGB2','B1:B11');
[alpha,sigma,GB2fit] = fn_recover_AcGB2(thetasp100gb2(1:4),thetasp100gb2(5:8),thetasp100gb2(10),thetasp100gb2(11),thetasp100gb2(9),Qsp100,nsp100);
for i=1:nt
GB2 =  generateGB2(1,1,thetasp100gb2(10),thetasp100gb2(11),nsp100);
[h,pGB2sp100(i)] = kstest2(GB2(burnin:nsp100),GB2fit(burnin:nsp100),'Alpha',0.01);
end


[
 sum(pFrsp100<0.01) sum(pFrsp100<0.05)  
  sum(pGB2sp100<0.01) sum(pGB2sp100<0.05)  
 
]
