function [gb2r] = generateGB2(a,b,p,q,n)

% rgb2(x, shape1, scale, shape2, shape3) GB2 package in R to be convenient.

U1 = chi2rnd(2*p,n,1);
U2 = chi2rnd(2*q,n,1);
Z = U1./U2;
Z = Z.^(1/a);
gb2r = b.*Z;

return
