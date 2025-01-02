function loglik = fn_loglik_AcGB2solver(theta)
% al 1x4 for sigma_t
% be 1x4 for alpha_t
n = evalin('base','n');
%p = evalin('base','p');
Qgb2 = evalin('base','Qgb2');
penalty = 10^10;
burnin = evalin('base','burnin');

mu = theta(9);
al = theta(1:4); %alpha0-3
be = theta(5:8); %beta0-3
p = theta(10);
q = theta(11);

la = zeros(n,1); lb = la;
% la(1) = al(1) - al(3)*exp(-mu*al(4));
% lb(1) = be(1) + be(3)*exp(-mu*be(4));
la(1) = (al(1)+al(3)/2)/(1-al(2));
lb(1) = (be(1)-be(3)/2)/(1-be(2));

for i=2:n
    la(i) = al(1) + al(2)*la(i-1) - al(3)*exp(-al(4)*Qgb2(i-1));
    lb(i) = be(1) + be(2)*lb(i-1) + be(3)*exp(-be(4)*Qgb2(i-1));
end
la = la(burnin+1:n);
lb = lb(burnin+1:n);
Qgb2lik = Qgb2(burnin+1:n);

a = exp(la);
b = exp(lb);
lk = -log(beta(p,q)) + mean(la - p*a.*lb + (p*a-1).*log(Qgb2lik-mu) - (p+q)*log(1+((Qgb2lik-mu)./b).^a));
if imag(lk)==0
    loglik = -lk; 
else
    loglik = penalty;
end
if loglik<-penalty
    loglik = penalty;
end

return

   