function y=mastereqn_optfun_lsq(p, data)
s=size(data);
sizet=s(3);
A=[0 1 0; 1 -4 1 ; 0 1 0];
%A=[1 2 1; 2 -12 2 ; 1 2 1];
y(1:(sizet-1))=0;
for t=1:(sizet-1)
    dcurrent=data(2,2,t);
    %dneighbour=sum(sum(data(:,:,t).*A));
    B=data(:,:,t).*A;
    dneighbour=sum(B(:));
    pred=dcurrent+p*dneighbour;
    y(t)=data(2,2,t+1)-pred;
end