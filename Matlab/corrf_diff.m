function g=corrf_diff(param,tau)
Amp=param(1);
tauD=param(2);
Offset=param(3);
g=Amp./(1+tau./tauD)+Offset;
