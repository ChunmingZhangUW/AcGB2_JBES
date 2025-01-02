function [a,b,GB2] = fn_recover_AcGB2(al,be,p,q,mu,Q,n)
% al 1x4 for sigma_t
% be 1x4 for alpha_t

 la = zeros(n,1); lb = la;
% la(1) = al(1) - al(3)*exp(-mu*al(4));
% lb(1) = be(1) + be(3)*exp(-mu*be(4));
 la(1) = (al(1)+al(3)/2)/(1-al(2));
 lb(1) = (be(1)-be(3)/2)/(1-be(2));
%la(1) = 0;
%lb(1) = 0;

for i=2:n
    la(i) = al(1) + al(2)*la(i-1) + al(3)*exp(-al(4)*Q(i-1));
    lb(i) = be(1) + be(2)*lb(i-1) - be(3)*exp(-be(4)*Q(i-1));
end

a = exp(la);
b = exp(lb);
GB2 = ((Q-mu)./b).^a;

return

   