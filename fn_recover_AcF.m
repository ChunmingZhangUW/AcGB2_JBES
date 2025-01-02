function [alpha,sigma,Fr] = fn_recover_AcF(gammap,betap,mu,Q,n)
% betap 1x4 for sigma_t
% gammap 1x4 for alpha_t

lsigma = zeros(n,1); lalpha = lsigma;
lalpha(1) = (gammap(1)+gammap(3)/2)/(1-gammap(2));
lsigma(1) = (betap(1)-betap(3)/2)/(1-betap(2));

for i=2:n
    lalpha(i) = gammap(1) + gammap(2)*lalpha(i-1) + gammap(3)*exp(-gammap(4)*Q(i-1));
    lsigma(i) = betap(1)  + betap(2)*lsigma(i-1)  -  betap(3)*exp(-betap(4)*Q(i-1)); 
end

alpha = exp(lalpha);
sigma = exp(lsigma);
Fr = ((Q-mu)./sigma).^alpha;
    
return

   