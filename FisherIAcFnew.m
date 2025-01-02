function [plog] = FisherIAcFnew(para,Q,burnin)

%#####################################################################
%### Function to Estimate the Fisher Information Matrix (FI) by using ###
%###   sample covariance matrix of 1st order partial derivatives   ###
%### Converted from R-cdoe made by Zifeng Zhao 05/29/2021, FisherFI revised 6/1/2021
%#####################################################################
  
  n=length(Q);
  gamma0 = para(2); gamma1 = para(3); gamma2 = para(4); gamma3 = para(5);
  beta0 = para(6); beta1 = para(7); beta2 = para(8); beta3 = para(9);
  mu = para(1);
    
  [alpha,sigma,loglik] = fn_loglik_AcF(para(2:5),para(6:9),mu,Q,n,burnin);
  k = burnin;
  et = 10^(-5); % tolorate error of gamma1^i  beta1^i
  
%   gb = max(gamma1, beta1);
%   k = ceil(log(et)/log(gb));
  
  gammaSeqSum = (1-gamma1^(k+1))/(1-gamma1);
  betaSeqSum = (1-beta1^(k+1))/(1-beta1);
  K = 0:k-1;
  gamma1K = gamma1.^K;
  beta1K = beta1.^K;
  Kr = k:-1:1;
  
   for i=k:n
      K1=((Q(i)-mu)/sigma(i))^(-alpha(i)); % equivalent to Y_i Frechet
      const1=1+log(K1)-K1*log(K1);  % page 7, \partial l_t(\theta)/\partial alpha_t = const1/alpha_t
      const2=alpha(i)*(1-K1);
      
      plog(1,i-k+1) = gammaSeqSum*const1;
      plog(2,i-k+1) = gamma1K*log(alpha(i-Kr+1))*const1;
      plog(3,i-k+1) = gamma1K*exp(-gamma3*Q(i-Kr+1))*const1;
      plog(4,i-k+1) = -gamma2*sum(beta1K.*Q(i-Kr+1)'.*exp(-gamma3*Q(i-Kr+1)'))*const1;
      
      plog(5,i-k+1) = betaSeqSum*const2;
      plog(6,i-k+1) = beta1K*log(sigma(i-Kr+1))*const2;
      plog(7,i-k+1) = beta1K*exp(-beta3*Q(i-Kr+1))*const2;
      plog(8,i-k+1) = beta2*sum(beta1K.*Q(i-Kr+1)'.*exp(-beta3*Q(i-Kr+1)'))*const2;
      
      plog(9,i-k+1)=(alpha(i)+1)/(Q(i)-mu)-alpha(i)/sigma(i).*((Q(i)-mu)/sigma(i))^(-alpha(i)-1);
   end
  
  
  return


