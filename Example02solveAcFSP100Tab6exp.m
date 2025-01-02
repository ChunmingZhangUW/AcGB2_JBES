
clear all;

rand('seed',5252921);

QFr = xlsread('SP100_processed16.csv', 'SP100_processed16', 'C2:C1258');
Para = xlsread('SP100paraacfNew.csv', 'SP100paraacfNew', 'B2:J401');
loglik = xlsread('SP100likeacfNew.csv', 'SP100likeacfNew', 'B2:B401');
n = length(QFr);
I = find(loglik<0);
loglik = -loglik(I);
Para = Para(I,:);
lb = [-10 0.01 0.001 1 -10 0.01 -10 1 0];
ub = [10 0.99 10 10 10 0.99 -0.01 10 0.981741135704124 ];
nI = length(I);
Is = [];
for i=1:nI
    pa = Para(i,:);
    pib = pa - lb;
    pub = pa - ub;
    I0l = find(pib==0);
    I0u = find(pub==0);
    if length(I0l)+length(I0u)==0
        Is = [Is i];
    end
end
loglik = loglik(Is);
Para = Para(Is,:);
[Mm I]= max(loglik);
 theta0 = Para(I,:);
 %return
theta0(7) = - theta0(7);
QFr = exp(QFr);
burnin = 100;
[alpha,sigma,loglik] = fn_loglik_AcF(theta0(1:4),theta0(5:8),theta0(9),QFr,n,burnin);

thetahat = [theta0(9) theta0(1:8)];
 [FInew] = FisherIAcFnew(thetahat,QFr,burnin);
  J=cov(FInew'); 
  
%#For computational stability
D12=diag(1./sqrt(diag(J)));
J1=D12*J*D12;
sdNormal=sqrt(diag(D12*inv(J1)*D12))/sqrt(n);
[theta0(9) sdNormal(9) 
    theta0(1:8)' sdNormal(1:8) 
    ]
loglik
%return
opfunc = @fn_loglik_AcFsolver;%

theta0 = [   -1.2529      0.74275       9.9448       1.6776     -0.33067      0.76588       4.0544 3.0329      0.91323]; % alpha_t, sigma_t, mu
[alpha0,sigma0,loglik0] = fn_loglik_AcF(theta0(1:4),theta0(5:8),theta0(9),QFr,n,burnin); %
lb = [ -10 0.01 0.001 1 -10 0.01 -10 1 0];
ub = [ 10 0.99 20 10 10 0.99 5 10 1];
loglik400=-inf; 
for i=1:10
theta1 = fmincon(opfunc,theta0,[],[],[],[],lb,ub);
[alpha1,sigma1,loglik1] = fn_loglik_AcF(theta1(1:4),theta1(5:8),theta1(9),QFr,n,burnin);
if loglik1>loglik400
    thetafinal = theta1;
    loglik400 = loglik1;
end
theta0 = thetafinal;
end
thetafinal
[alpha400,sigma400,loglik401] = fn_loglik_AcF(thetafinal(1:4),thetafinal(5:8),thetafinal(9),QFr,n,burnin);

thetahat = [thetafinal(9) thetafinal(1:8)];
 [FInew] = FisherIAcFnew(thetahat,QFr,burnin);
  J=cov(FInew'); 
  
%#For computational stability
D12=diag(1./sqrt(diag(J)));
J1=D12*J*D12;
sdNormal=sqrt(diag(D12*inv(J1)*D12))/sqrt(n);
[thetafinal(9) sdNormal(9) 
    thetafinal(1:8)' sdNormal(1:8) 
    ]
loglik401