library(magrittr)
library(ggplot2)
nfdjia<-read.csv("djia.csv")
nfdjia$Date<-as.Date(nfdjia$Date)
stockrld<-nfdjia$Maxima
stockrld1<-nfdjia


stockrld1%>%
  ggplot(aes(x=Date,y=Maxima))+
  geom_point()+
  geom_smooth()
