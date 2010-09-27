function y=mastereqn_optfun(p, data, dt)
y=0;
s=size(data);
sizet=s(3);
A=[0 1 0; 1 -4 1 ; 0 1 0];
for t=1:(sizet-1)
    dcurrent=data(2,2,t);
    B=data(:,:,t).*A;
    dneighbour=sum(B(:));
    pred=dcurrent+p*dt*dneighbour;
    err=(data(2,2,t+1)-pred)^2;
    y=y+err;
end