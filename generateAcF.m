function [Q,a,b] = generateAcF(theta,n)

% a_t - nx1
% b_t - nx1
% theta = (alphap(1-4), betap(1-4),mu)

mu = theta(9);
al = theta(1:4); %gamma0-3
be = theta(5:8); %beta0-3
burnin = 10000;
Fr = 1./exprnd(1,n+burnin,1);

la = zeros(n+burnin,1); lb = la; a = la; b = a; Q = a;
% la(1) = al(1) - al(3)*exp(-mu*al(4));
% lb(1) = be(1) + be(3)*exp(-mu*be(4));
lal(1) = (al(1)+al(3)/2)/(1-al(2));
lsi(1) = (be(1)-be(3)/2)/(1-be(2));

Q(1) = 1;
for i=2:n+burnin
    lal(i) = al(1) + al(2)*la(i-1) + al(3)*exp(-al(4)*Q(i-1));  % log(alpha_t) in the paper
    lsi(i) = be(1) + be(2)*lb(i-1) - be(3)*exp(-be(4)*Q(i-1));  % log(sigma_t) in the paper
    a(i) = exp(lal(i));
    b(i) = exp(lsi(i));
    Q(i) = mu + b(i)*Fr(i)^(1/a(i));
end

a = a(burnin+1:n+burnin);
b = b(burnin+1:n+burnin);
Q = Q(burnin+1:n+burnin);

return

