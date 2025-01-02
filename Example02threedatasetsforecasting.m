
clear all;
clf

%rng('default')
rand('seed',5252921);
burnin = 100; 
ns = 1000; %simulation sample size
ni = 0; %recovered sample size 
n975 = floor((ni+ns)*0.95);
n25 = floor((ni+ns)*0.05);

Qsp100 = xlsread('SP100_processed.csv', 'SP100_processed', 'C2:C1258');
nsp100 = length(Qsp100); Qsp100 = exp(Qsp100);
Qsp100pr = xlsread('sp100_prediction.csv', 'sp100_prediction', 'C2:C736');
nsp100pr = length(Qsp100pr); Qsp100pr = exp(Qsp100pr);
[~,xdates,~] = xlsread('sp100_prediction.csv', 'sp100_prediction', 'B2:B736');
date=datestr(xdates,'mm/dd/yyyy');

thetasp100acf = xlsread('ThreeDatasetsparameterestimations.xlsx','AcF','F1:F9');
[alpha,sigma,Frfit] = fn_recover_AcF(thetasp100acf(1:4),thetasp100acf(5:8),thetasp100acf(9),Qsp100,nsp100);
[Q25, Q975] = AcFforecast(nsp100,alpha,sigma,ni,Frfit, Qsp100,Qsp100pr,nsp100pr,ns,thetasp100acf,n975,n25);
figure(2)
hold on
plot(Q25,'m')
I = 1:40:nsp100pr;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 735.5]);
plot(Qsp100pr,'b')
plot(Q975,'r')
hold off
ylabel('AcF SP100 forecast','fontsize',12)
ylim([1,2.3])
[sum(Qsp100pr>Q975) sum(Qsp100pr>Q975)/nsp100pr max(Q975)]
[sum(Qsp100pr>Q25) 1-sum(Qsp100pr>Q25)/nsp100pr]
%return


thetasp100gb2 = xlsread('ThreeDatasetsparameterestimations.xlsx','AcGB2','B1:B11');
[alpha,sigma,GB2fit] = fn_recover_AcGB2(thetasp100gb2(1:4),thetasp100gb2(5:8),thetasp100gb2(10),thetasp100gb2(11),thetasp100gb2(9),Qsp100,nsp100);
[Q25, Q975] = AcGB2forecast(nsp100,alpha,sigma,ni,GB2fit, Qsp100,Qsp100pr,nsp100pr,ns,thetasp100gb2,n975,n25);
figure(5)
hold on
plot(Q25,'m')
I = 1:40:nsp100pr;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 735.5]);
plot(Qsp100pr,'b')
plot(Q975,'r')
hold off
ylabel('AcGB2 SP100 forecast','fontsize',12)
ylim([1,2.3])

[sum(Qsp100pr>Q975) sum(Qsp100pr>Q975)/nsp100pr max(Q975)]
[sum(Qsp100pr>Q25) 1-sum(Qsp100pr>Q25)/nsp100pr]
%return


return