function [Q,a,b] = generateAcGB2(theta,n)


mu = theta(9);
al = theta(1:4); 
be = theta(5:8); 
p = theta(10);
q = theta(11);
burnin = 10000;
gb2r = generateGB2(1,1,p,q,n+burnin);
%gb2r(1:10)
la = zeros(n+burnin,1); lb = la; a = la; b = a; Q = a;
% la(1) = al(1) - al(3)*exp(-mu*al(4));
% lb(1) = be(1) + be(3)*exp(-mu*be(4));
la(1) = (al(1)+al(3)/2)/(1-al(2));
lb(1) = (be(1)-be(3)/2)/(1-be(2));

Q(1) = 1; %1 before
for i=2:n+burnin
    la(i) = al(1) + al(2)*la(i-1) - al(3)*exp(-al(4)*Q(i-1));  % sigma_t in the paper
    lb(i) = be(1) + be(2)*lb(i-1) + be(3)*exp(-be(4)*Q(i-1));  % alpha_t in the paper
    a(i) = exp(la(i));
    b(i) = exp(lb(i));
    Q(i) = mu + a(i)*gb2r(i)^(1/b(i));
end

a = a(burnin+1:n+burnin);
b = b(burnin+1:n+burnin);
Q = Q(burnin+1:n+burnin);
%Q(1:10)
return

