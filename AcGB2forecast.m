function [Q25, Q975] = AcGB2forecast(ndjia,alpha,sigma,ni,GBfit, Qdjia,Qdjiapr,ndjiapr,ns,thetadjiaGB2,n975,n25)
a975 = alpha(ndjia); b975=sigma(ndjia);
a25 = alpha(ndjia); b25=sigma(ndjia);
GBi975 =[]; GBi25 = GBi975;
if ni>0
GBi975 = GBfit(ndjia-ni+1:ndjia); GBi25 = GBi975;
end
Qtm1 = Qdjia(ndjia);
for i=1:ndjiapr
GBs = generateGB2(1,1,thetadjiaGB2(10),thetadjiaGB2(11),ns); 
GB975 =[GBi975' GBs']';
GB975 = sort(GB975);
GB975t = GB975(n975);
[Q975(i,1),a975,b975] = forecastAcGB2(thetadjiaGB2, log(a975), log(b975), Qtm1, GB975t);
if ni>0
GBi975 = [GBi975(2:ni)' GB975t]';
else
    GBi975 =[];
end
GB25 =[GBi25' GBs']';
GB25 = sort(GB25);
GB25t = GB25(n25);
[Q25(i,1),a25,b25] = forecastAcGB2(thetadjiaGB2, log(a25), log(b25), Qtm1, GB25t);
if ni>0
GBi25 = [GBi25(2:ni)' GB25t]';
else
    GBi25 = [];
end
Qtm1 = Qdjiapr(i);
end