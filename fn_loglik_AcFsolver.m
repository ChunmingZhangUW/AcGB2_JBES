function loglik = fn_loglik_AcFsolver(theta)
% al 1x4 for sigma_t
% be 1x4 for alpha_t
n = evalin('base','n');
%p = evalin('base','p');
QFr = evalin('base','QFr');
burnin = evalin('base','burnin');

mu = theta(9);
al = theta(1:4); % alpha_t
be = theta(5:8); % sigma_t

la = zeros(n,1); lb = la;
% la(1) = al(1) - al(3)*exp(-mu*al(4));
% lb(1) = be(1) + be(3)*exp(-mu*be(4));
la(1) = (al(1)+al(3)/2)/(1-al(2));
lb(1) = (be(1)-be(3)/2)/(1-be(2));


for i=2:n
    la(i) = al(1) + al(2)*la(i-1) + al(3)*exp(-al(4)*QFr(i-1)); %log_alpha_t
    lb(i) = be(1) + be(2)*lb(i-1) - be(3)*exp(-be(4)*QFr(i-1)); %log_sigma_t
    %[i la(i) lb(i) QFr(i-1)]
end
la = la(burnin+1:n);
lb = lb(burnin+1:n);
QFrlik = QFr(burnin+1:n);

a = exp(la);
b = exp(lb);
lk = mean(la + a.*lb - (a+1).*log(QFrlik-mu) - b.^a.*(QFrlik-mu).^(-a));
if imag(lk)==0
    loglik = -lk; % the smaller the better, different from loglik
else
    loglik = 10^10;
end

return

   