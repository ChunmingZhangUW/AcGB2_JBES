function [Q,a,b] = forecastAcGB2(theta, la, lb, Qtm1, GB2)

mu = theta(9);
al = theta(1:4); %gamma0-3
be = theta(5:8); %beta0-3
p = theta(10);
q = theta(11);
%gb2r = generateGB2(1,1,p,q,1);

    lal = al(1) + al(2)*la + al(3)*exp(-al(4)*Qtm1);  % log(alpha_t) in the paper
    lsi = be(1) + be(2)*lb - be(3)*exp(-be(4)*Qtm1);  % log(sigma_t) in the paper
    a = exp(lal);
    b = exp(lsi);
    Q = mu + b*GB2^(1/a);


return

