
clear all;

rand('seed',5252921);

Qgb2 = xlsread('SP100_processed.csv', 'SP100_processed', 'C2:C1258');
Para = xlsread('SP100paraacgb2New.csv', 'SP100paraacgb2New', 'B2:L401');
loglik = xlsread('SP100likeacgb2New.csv', 'SP100likeacgb2New', 'B2:B401');
n = length(Qgb2);
I = find(loglik<0);
loglik = -loglik(I);
Para = Para(I,:);
lb = [-0.5 0.70 17.0 2.5 -0.50 0.6 -6.0 4.0 0.70 1 0.1];
ub = [-0.4 0.90 18.0 3.0 -0.20 0.8 -5.0 6.0 0.95 2 0.5];
nI = length(I);
Is = [];
for i=1:nI
    pa = Para(i,:);
    pib = pa - lb;
    pub = pa - ub;
    I0l = find(pib==0);
    I0u = find(pub==0);
    if length(I0l)+length(I0u)<=5
        Is = [Is i];
    end
end
loglik = loglik(Is);
Para = Para(Is,:);
[Mm I]= max(loglik);
 theta0 = Para(I,:);
 %return
theta0(7) = - theta0(7);
Qgb2 = exp(Qgb2);
burnin = 100;
[alpha,sigma,loglik] = fn_loglik_AcGB2(theta0(1:4),theta0(5:8),theta0(10), theta0(11), theta0(9),Qgb2,n,burnin);

thetahat = [theta0(9) theta0(1:8) theta0(10:11)];
 [FInew] = FisherIAcGB2new(thetahat,Qgb2,burnin);
  J=cov(FInew'); 
  
%#For computational stability
D12=diag(1./sqrt(diag(J)));
J1=D12*J*D12;
sdNormal=sqrt(diag(D12*inv(J1)*D12))/sqrt(n);
[theta0' sdNormal]
loglik
return
