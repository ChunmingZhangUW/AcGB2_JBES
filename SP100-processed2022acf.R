library(tidyverse)
library(stringr)
library(GB2)
library(extraDistr)
library(stats)
nfsp100<-read.csv("SP100_processed.csv")
nfsp100$Date<-as.Date(nfsp100$Date)
stockrlddd<-nfsp100$Maxima
stockrlddd1<-nfsp100

R<-exp(stockrlddd)
set.seed(632021)
#use acgb2 to est sp100


llk6 <- function(theta, R, ub, lb, burnin=100){
  
  n=length(R)
  alat=theta[10]
  bett=theta[11]
  pt=rep(0,n) #To store estimated alpha_t  
  qt=rep(0,n) #To store estimated sigma_t  
  
  #Recover latent alpha and sigma
  pt <- rec_pt(theta, R)
  qt <- rec_qt(theta, R)
  
  
  ###########################################
  #penalize the theta not in parameter space#
  ###########################################
  penal=10^10;penalty=penal*(sum(theta<lb)+sum(theta>ub)) 
  
  #Truncate the first burnin of the data for better estimation of alpha_t
  # R1=(R-theta[9])/est_sigma
  pt=pt[burnin:n]
  qt=qt[burnin:n]
  R=R[burnin:n]
  # loglik=-1/n*sum(log(pt)-log(qt)-(pt+1)*log((R)/qt)-((R)/qt)^(3*pt))
  loglik=-1/n*sum(log(pt)-log(beta(alat,bett))+(pt*alat-1)*log(R-theta[9])-(pt*alat)*log(qt)-(alat+bett)*log(1+((R-theta[9])/qt)^pt))
  # loglik=-1/n*sum()
  # llkin[i]<-log(pt[i])-log(beta(alat[i],bett[i]))+(pt[i]*alat[i]-1)*log(gbt[i])-(pt[i]*alat[i])*log(qt[i])-(alat[i]+bett[i])*log(1+(gbt[i]/qt[i]))
  return(loglik+penalty)
}

rec_pt<-function(para,gbt){
  
  b0=para[1];b1=para[2];b2=para[3];b3=para[4]
  n=length(gbt)
  pt=rep(0,n)
  pt[1]=exp((b0+b2/2)/(1-b1))
  for(i in 1:(n-1)){
    pt[i+1]=exp(b0+b1*log(pt[i])+b2*exp(-b3*gbt[i]))
  }
  return(pt)
}

rec_qt<-function(para,gbt){
  
  r0=para[5];r1=para[6];r2=para[7];r3=para[8]
  n=length(gbt)
  qt=rep(0,n)
  qt[1]=exp((r0+r2/2)/(1-r1))
  for(i in 1:(n-1)){
    qt[i+1]=exp(r0+r1*log(qt[i])+r2*exp(-r3*gbt[i]))
  }
  return(qt)
}

resultm<-matrix(rep(0,(11*400)),nrow=400,ncol=11)
likeacgb2<-numeric(400)


for (i in 1:400)
{
  lb1=c(-1,0.01,0.01,1,-2,0.01,-1,1,0,0,0)
  ub1=c(1, 0.99,1,   10,2,0.99,-0.01,10,0.981741135704124,10,10)
  loc1 <- runif(1,0,1)
  start1 <- loc1*lb1 + (1-loc1)*ub1
  fitTMP=nlminb(start1,llk6,lower=lb1,upper=ub1,control=list(eval.max=1e7, iter.max=1e7),
                R=R,ub=ub1,lb=lb1,burnin=100)
  likeacgb2[i]=fitTMP$objective
  resultm[i,]=fitTMP$par
}
#find mean of the 50 sets of parameters
resultmd<-data.frame(resultm)

resultmd<-resultmd%>%
  filter(X3!=1&X6!=0.01&X2!=X6&X2>0.8)
resultf<-rep(0,11)
for (i in 1:11)
{
  resultf[i]<-mean(resultmd[,i])
}

signif(resultf,3)

vresultf<-rep(0,11)

for (i in 1:11)
{
  vresultf[i]<-sd(resultmd[,i])
}

signif(vresultf,3)






#use acf to est sp100


rec_alphaLog<-function(para,R){
  omega=para[1];b=para[2];c=para[3];d=para[4]
  n=length(R)
  
  logAlpha=rep(0,n)
  logAlpha[1]=(omega+c/2)/(1-b)
  
  #recover the latent alpha series
  for(i in 2:n){
    logAlpha[i]=omega+b*logAlpha[i-1]+c*exp(-d*R[i-1])
  }
  return(exp(logAlpha))
}

rec_sigma<-function(para,R){
  
  beta0=para[5];beta1=para[6];beta2=para[7];beta3=para[8]
  n=length(R)
  
  logSig=rep(0,n)
  logSig[1]=(beta0+beta2/2)/(1-beta1)
  
  #recover the latent alpha series
  for(i in 2:n){
    logSig[i]=beta0+beta1*logSig[i-1]+beta2*exp(-beta3*R[i-1])
  }
  return(exp(logSig))
}

fn_loglik <- function(theta, R, ub, lb, burnin=100){
  
  n=length(R)
  est_alpha=rep(0,n) #To store estimated alpha_t  
  est_sigma=rep(0,n) #To store estimated sigma_t  
  
  #Recover latent alpha and sigma
  est_alpha <- rec_alphaLog(theta, R)
  est_sigma <- rec_sigma(theta, R)
  
  ###########################################
  #penalize the theta not in parameter space#
  ###########################################
  penal=10^10;penalty=penal*(sum(theta<lb)+sum(theta>ub)) 
  
  #Truncate the first burnin of the data for better estimation of alpha_t
  R1=(R-theta[9])/est_sigma
  est_alpha=est_alpha[burnin:n]
  est_sigma=est_sigma[burnin:n]
  R1=R1[burnin:n]
  
  #Loglikelihood for GEV model
  loglik=-1/n*sum(log(est_alpha)-log(est_sigma)-(est_alpha+1)*log(R1)-(R1)^(-est_alpha))
  return(loglik+penalty)      
}





resultm2<-matrix(rep(0,(9*400)),nrow=400,ncol=9)
likeacf<-numeric(400)

for (i in 1:400)
{
  lb1=c(-2,0.01,0.001, 1,-2,0.01,  -10,1, 0)
  ub1=c( 2,0.99,   10,10, 2,0.99,-0.01,10,0.981741135704124)
  loc1 <- runif(1,0,1)
  start1 <- loc1*lb1 + (1-loc1)*ub1
  fitTMP=nlminb(start1,fn_loglik,lower=lb1,upper=ub1,control=list(eval.max=1e7, iter.max=1e7),
                R=R,ub=ub1,lb=lb1,burnin=100)
  likeacf[i]=fitTMP$objective
  resultm2[i,]=fitTMP$par
}



resultmd2<-data.frame(resultm2)

#resultmd2<-resultmd2%>%
#  dplyr::filter(X2!=0.01&X4!=1&X6!=0.01&X7!=-0.01&X8!=1&X1!=X5&X1<0&X1>(-0.3)&X8<8)
#resultmd2<-resultmd2%>%
#  dplyr::filter(X2!=0.01&X4!=1&X6!=0.01&X7!=-0.01&X8!=1&X1!=X5&X4<10&X6>0.03&X6<0.04)

resultf2<-rep(0,9)

for (i in 1:9)
{
  resultf2[i]<-mean(resultmd2[,i])
}
signif(resultf2,3)

vresultf2<-rep(0,9)

for (i in 1:9)
{
  vresultf2[i]<-sd(resultmd2[,i])
}

signif(vresultf2,3)

# write.csv(resultm, file = !"Tab6paraacgb2New.csv")
# write.csv(likeacgb2, file =! "Tab6likeacgb2New.csv")
# write.csv(resultm2, file = "SP100paraacfNew.csv")
# write.csv(likeacf, file = "SP100likeacfNew.csv")
