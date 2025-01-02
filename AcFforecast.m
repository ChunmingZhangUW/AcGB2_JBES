function [Q25, Q975] = AcFforecast(ndjia,alpha,sigma,ni,Frfit, Qdjia,Qdjiapr,ndjiapr,ns,thetadjiaacf,n975,n25)
a975 = alpha(ndjia); b975=sigma(ndjia);
a25 = alpha(ndjia); b25=sigma(ndjia);
Fri975 =[]; Fri25 = Fri975;
if ni>0
Fri975 = Frfit(ndjia-ni+1:ndjia); Fri25 = Fri975;
end
Qtm1 = Qdjia(ndjia);
for i=1:ndjiapr
Frs = 1./exprnd(1,ns,1); 
Fr975 =[Fri975' Frs']';
Fr975 = sort(Fr975);
Fr975t = Fr975(n975);
[Q975(i,1),a975,b975] = forecastAcF(thetadjiaacf, log(a975), log(b975), Qtm1, Fr975t);
if ni>0
Fri975 = [Fri975(2:ni)' Fr975t]';
else
    Fri975 =[];
end
Fr25 =[Fri25' Frs']';
Fr25 = sort(Fr25);
Fr25t = Fr25(n25);
[Q25(i,1),a25,b25] = forecastAcF(thetadjiaacf, log(a25), log(b25), Qtm1, Fr25t);
if ni>0
Fri25 = [Fri25(2:ni)' Fr25t]';
else
    Fri25 = [];
end
Qtm1 = Qdjiapr(i);
end