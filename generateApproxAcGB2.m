function [Q,a,b] = generateApproxAcGB2(theta,Sigma2,n,m)

% m(1) - dimension
% m(2) - GB2 dimension
% m(3) - t4 dimension
% m(4) - normal dimension

mu = theta(9);
al = theta(1:4); 
be = theta(5:8); 
p = theta(10);
q = theta(11);
rng('default')  % For reproducibility

mu0 = zeros(m(1),1);
R = mvnrnd(mu0,Sigma2,n+100);
%p = mvncdf(R, mu0, Sigma2);
%gb2r = generateGB2(1,1,p,q,n+100);
a=1; b=1;
gb2ri = zeros(n+100,m(1));

for i=1:m(2)
%U1 = chi2rnd(2*p,n,1);
P = normcdf(R(:,i),0,sqrt(Sigma2(i,i)));
U1 = chi2inv(P,2*p);
U2 = chi2rnd(2*q,n+100,1);
Z = U1./U2;
Z = Z.^(1/a);
gb2ri(:,i) = b.*Z;
end

for i=m(2)+1:m(2)+m(3)
    P = normcdf(R(:,i),0,sqrt(Sigma2(i,i)));
    gb2ri(:,i) = exp(tinv(P,6));
end

for i=m(2)+m(3)+1:sum(m(2:4))
    gb2ri(:,i) = exp(R(:,i));
end
gb2r = max(gb2ri,[],2);
%return
la = zeros(n+100,1); lb = la; a = la; b = a; Q = a;
% la(1) = al(1) - al(3)*exp(-mu*al(4));
% lb(1) = be(1) + be(3)*exp(-mu*be(4));
la(1) = (al(1)+al(3)/2)/(1-al(2));
lb(1) = (be(1)-be(3)/2)/(1-be(2));

Q(1) = 1;
for i=2:n+100
    la(i) = al(1) + al(2)*la(i-1) - al(3)*exp(-al(4)*Q(i-1));  % sigma_t in the paper
    lb(i) = be(1) + be(2)*lb(i-1) + be(3)*exp(-be(4)*Q(i-1));  % alpha_t in the paper
    a(i) = exp(la(i));
    b(i) = exp(lb(i));
    Q(i) = mu + a(i)*gb2r(i)^(1/b(i));
end

a = a(101:n+100);
b = b(101:n+100);
Q = Q(101:n+100);

return

