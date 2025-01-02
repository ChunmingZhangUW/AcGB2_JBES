clear all;

Qsp100 = xlsread('sp100_processed.csv', 'sp100_processed', 'C2:C1258');
nsp100 = length(Qsp100); Qsp100 = exp(Qsp100);
[~,xdates,~] = xlsread('sp100_processed.csv', 'sp100_processed', 'B2:B1258');
date=datestr(xdates,'mm/dd/yyyy');

figure(1)
plot(Qsp100,'k.')
hold on
I = 1:60:nsp100;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 1257.5]);
hold off
ylabel('Exp(maxima)','fontsize',12)


[sum(Qsp100<1)/nsp100 sum(Qsp100>exp(0.05))/nsp100 sum(Qsp100>exp(0.1))/nsp100 sum(Qsp100>exp(0.15))/nsp100 sum(Qsp100>exp(0.20))/nsp100 sum(Qsp100>exp(0.3))/nsp100]
A=[min(Qsp100) mean(Qsp100) median(Qsp100) max(Qsp100)]
log(A)
