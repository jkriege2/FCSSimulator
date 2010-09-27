function y=mastereqn_optfun1_lsq(p, data, dt)
s=size(data);
sizet=s(3);
A1=[0 1 0; 0 -1 0 ; 0 0 0];
A2=[0 0 0; 1 -1 0 ; 0 0 0];
A3=[0 0 0; 0 -1 1 ; 0 0 0];
A4=[0 0 0; 0 -1 0 ; 0 1 0];
y(1:ceil((sizet-1)/4))=0;
for t=4:4:(sizet-1)
    B=data(:,:,t).*A1+data(:,:,t-1).*A2+data(:,:,t-2).*A3+data(:,:,t-3).*A4;
    dneighbour=p*dt*sum(B(:));
    pred=data(2,2,t)+dneighbour;
    y(t/4)=data(2,2,t+1)-pred;
end