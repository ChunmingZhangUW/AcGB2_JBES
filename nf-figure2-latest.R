# CDF and PDF

####################################
# Different number of time series for AcF
####################################

library(GB2)
library(fitdistrplus)
library(evd)

max1000<-numeric(10000)
for(i in 1:10000){
  max1000[i]<-max(qgb2(runif(1000,0,1)^1000,2,5,3,3))}
#plot(density(max1000),xlim=c(0,50),ylim=c(0,0.3))
#x<-seq(0,50,length=1000)
#y<-dgb2(x,2,5,3,3)
#lines(x,y,col='red')

max900<-numeric(10000)
for(i in 1:10000){
  max900[i]<-max(qgb2(runif(900,0,1)^1000,2,5,3,3))}

max700<-numeric(10000)
for(i in 1:10000){
  max700[i]<-max(qgb2(runif(700,0,1)^1000,2,5,3,3))}


max500<--numeric(10000)
for(i in 1:10000){
  max500[i]<-max(qgb2(runif(500,0,1)^1000,2,5,3,3))}

fitfrechet1000<-fitdist(max1000, "frechet", method = "mle", lower = c(0, 0, 0), start = list(loc=0.1,scale=1, shape=1))


fitfrechet900<-fitdist(max900, "frechet", method = "mle", lower = c(0, 0, 0), start = list(loc=0.1,scale=1, shape=1))
fitfrechet700<-fitdist(max700, "frechet", method = "mle", lower = c(0, 0, 0), start = list(loc=0.1,scale=1, shape=1))
fitfrechet500<-fitdist(max500, "frechet", method = "mle", lower = c(0, 0, 0), start = list(loc=0.1,scale=1, shape=1))
#print(fitfrechet1000)
#print(fitfrechet900)
#print(fitfrechet700)
#print(fitfrechet500)




par(mfrow=c(1,2))
set.seed(0)

x_value=seq(0,30,length=1000)
x_value
max900f<-dfrechet(x_value,3.175270e-08, 3.809034e+00 ,2.128417e+00 )
max700f<-dfrechet(x_value,5.396968e-08 ,3.347004e+00  ,1.886508e+00    )
max500f<-dfrechet(x_value,7.760515e-09 , 2.728293e+00  ,1.493496e+00   )

plot(x_value,max900f,xlim=c(0,30),ylim=c(0,0.25),type='l',lty=4,xlab='',ylab='',main='')
lines(x_value,max700f,lty=2)
lines(x_value,max500f,lty=3)
y<-dgb2(x_value,2,5,3,3)
lines(x_value,y,lty=1)
legend("topright",legend=c("m=500","m=700","m=900","GB2(2,5,3,3)"),lty=c(3,2,4,1),cex=0.6)

x_value<-seq(0,30,length=1000)
y<-pgb2(x,2,5,3,3)
max900f<-pfrechet(x_value,3.175270e-08, 3.809034e+00 ,2.128417e+00 )
max700f<-pfrechet(x_value,5.396968e-08 ,3.347004e+00  ,1.886508e+00    )
max500f<-pfrechet(x_value,7.760515e-09 , 2.728293e+00  ,1.493496e+00   )

plot(x_value,max900f,ylim=c(0,1.1),type='l',lty=4,xlab='',ylab='',main='')
lines(x_value,max700f,lty=2)
lines(x_value,max500f,lty=3)
lines(x_value,y,lty=1)
legend("bottomright",legend=c("m=500","m=700","m=900","GB2(2,5,3,3)"),lty=c(3,2,4,1),cex=0.6)
