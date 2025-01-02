clear all;

Qsp100 = xlsread('sp100_processed.csv', 'sp100_processed', 'D101:I1258');
[~,xdates,~] = xlsread('sp100_processed.csv', 'sp100_processed', 'B101:B1258');
date=datestr(xdates,'mm/dd/yyyy');
nsp100 = length(Qsp100(:,1));
meanV = xlsread('sp100_combine_v1.csv', 'sp100_combine_v1', 'CZ101:CZ1258');
meanV = (meanV - mean(meanV))/std(meanV);
Qsp100(:,1:2)=Qsp100(:,5:6);


figure(1)
plot(Qsp100(:,1),'k.')
hold on
I = 1:60:nsp100;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 1157.5]);
hold off
ylabel('AcF \alpha_t','fontsize',12)

figure(2)
plot(Qsp100(:,2),'k.')
hold on
I = 1:60:nsp100;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 1157.5]);
hold off
ylabel('AcF \sigma_t','fontsize',12)

figure(3)
plot(Qsp100(:,3),'k.')
hold on
I = 1:60:nsp100;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 1157.5]);
hold off
ylabel('AcGB2 A_t','fontsize',12)

figure(4)
plot(Qsp100(:,4),'k.')
hold on
I = 1:60:nsp100;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 1157.5]);
hold off
ylabel('AcGB2 B_t','fontsize',12)

figure(5)
AcF = 1./Qsp100(:,1);
AcF = (AcF- mean(AcF))/std(AcF);
plot(AcF,'k')
hold on
plot(meanV,'r')
I = 1:60:nsp100;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 1157.5]);
hold off
ylabel('AcF vs GARCH','fontsize',12)
corr(AcF,meanV)

figure(6)
ACGB2 = 1./Qsp100(:,3);
AcGB2 = (ACGB2- mean(ACGB2))/std(ACGB2);
plot(AcGB2,'k')
hold on
plot(meanV,'r')
I = 1:60:nsp100;
set(gca,'xtick',I);
set(gca,'xticklabel',date(I,:),'fontsize',8)
xlim([0.8, 1157.5]);
hold off
ylabel('AcGB2 vs GARCH','fontsize',12)
corr(AcGB2,meanV)