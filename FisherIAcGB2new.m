function [plog] = FisherIAcGB2(para,Q,burnin)

%#####################################################################
%### Function to Estimate the Fisher Information Matrix (FI) by using ###
%###   sample covariance matrix of 1st order partial derivatives   ###
%### Converted from R-cdoe made by Zifeng Zhao 05/29/2021
%#####################################################################
  
% R = QFr

  n=length(Q);
  gamma0 = para(2); gamma1 = para(3); gamma2 = para(4); gamma3 = para(5);
  beta0 = para(6); beta1 = para(7); beta2 = para(8); beta3 = para(9);
  p = para(10); q = para(11);
  mu = para(1);
  
  %#Recover alpha_t from the data
  %alpha=rec_alphaLog(para,R);
  %sigma=rec_sigma(para,R);
  [alpha,sigma,loglik] = fn_loglik_AcGB2(para(2:5),para(6:9),p,q,mu,Q,n,burnin);
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
  
  
  %#Calculate the partial derivative for 9 parameters
  for i=k:n
      %[Q(i) mu sigma(i) alpha(i)]
      K1=((Q(i)-mu)/sigma(i))^(-alpha(i)); % equivalent to GB2 
      const1=1+p*log(K1)-(p+q)*(K1/(1+K1))*log(K1);  % 
      const2=alpha(i)*((p+q)*K1/(1+K1)-p);
      if imag(K1)~=0
          [Q(i) mu sigma(i) alpha(i)]
      pause
      end
      plog(1,i-k+1) = gammaSeqSum*const1;
      plog(2,i-k+1) = gamma1K*log(alpha(i-Kr+1))*const1;
      plog(3,i-k+1) = gamma1K*exp(-gamma3*Q(i-Kr+1))*const1;
      plog(4,i-k+1) = -gamma2*sum(beta1K.*Q(i-Kr+1)'.*exp(-gamma3*Q(i-Kr+1)'))*const1;
      
      plog(5,i-k+1) = betaSeqSum*const2;
      plog(6,i-k+1) = beta1K*log(sigma(i-Kr+1))*const2;
      plog(7,i-k+1) = beta1K*exp(-beta3*Q(i-Kr+1))*const2;
      plog(8,i-k+1) = beta2*sum(beta1K.*Q(i-Kr+1)'.*exp(-beta3*Q(i-Kr+1)'))*const2;
      
      %#Partial derivative for mu
      plog(9,i-k+1) = ((p+q)*alpha(i)*K1/(1+K1)-(alpha(i)*p-1))/(Q(i)-mu);
      plog(10,i-k+1)=-(psi(p)-psi(p+q))+log(K1)-log(1+K1);
      plog(11,i-k+1)=-(psi(q)-psi(p+q))-log(1+K1);
  end
  
  return


