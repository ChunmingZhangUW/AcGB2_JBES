function [alpha,sigma,loglik] = fn_loglik_AcF(gammap,betap,mu,Q,n,burnin)
% betap 1x4 for sigma_t
% gammap 1x4 for alpha_t
% if nargin<8
%     burnin = 1;
% end
 I = burnin:1:n;

lsigma = zeros(n,1); lalpha = lsigma;
lalpha(1) = (gammap(1)+gammap(3)/2)/(1-gammap(2)); %output NaN for DJI30
lsigma(1) = (betap(1)-betap(3)/2)/(1-betap(2));
 %lalpha(1) = 0;
 %lsigma(1) = 0;

for i=2:n
    lalpha(i) = gammap(1) + gammap(2)*lalpha(i-1) + gammap(3)*exp(-gammap(4)*Q(i-1));
    lsigma(i) = betap(1)  + betap(2)*lsigma(i-1)  -  betap(3)*exp(-betap(4)*Q(i-1)); 
end

alpha = exp(lalpha); %alpha(1:10) matched with R
sigma = exp(lsigma); % sigma(1:10) matched
loglik = mean(lalpha(I) + alpha(I).*lsigma(I) - (alpha(I)+1).*log(Q(I)-mu) - sigma(I).^alpha(I).*(Q(I)-mu).^(-alpha(I)));
return

   