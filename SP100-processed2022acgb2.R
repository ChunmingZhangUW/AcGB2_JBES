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
  #lb1=c(-0.9821,0.8374,8.7166,1.8524,-0.2911,0.7621,-4.2230,2.6699,0.9314,0,0)
  #ub1=c(-0.9819,0.8376,8.7168,1.8526,-0.2909,0.7623,-4.2228,2.6701,0.9316,10,10)
#  lb1=c(-1.2,0.70,7.0,1.6,-0.40,0.6,-5.0,2.4,0.83,0,0)
#  ub1=c(-0.8,0.90,9.0,2.0,-0.15,0.9,-4.0,2.9,0.95,10,10)
  #lb1=c(-0.5,0.70,17.0, 2.5,-0.50,0.6,-6.0,4.0,0.70,1,0.1)
  #ub1=c(-0.4,0.90,18.0, 3.0,-0.20,0.8,-5.0,6.0,0.95,2,0.5)
  lb1=c(-2.2,0.60,10,1.5, -1.50,0.10,-6.0, 3.0,0.7,0.7,0.1)
  ub1=c(-0.2,0.90,15,2.0, -1.00,0.60,-4.0, 5.0,0.9,0.95,0.3)
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




# write.csv(resultm, file = "SP100paraacgb2New.csv")
# write.csv(likeacgb2, file = "SP100likeacgb2New.csv")
# write.csv(resultm2, file =! "Tab6paraacfNew.csv")
# write.csv(likeacf, file = !"Tab6likeacfNew.csv")
