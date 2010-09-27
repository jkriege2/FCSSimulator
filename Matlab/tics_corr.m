function [corrf] = tics_corr(data)

ds=size(data);
N=ds(3);
XS=ds(1);
YS=ds(2);
corrf=1:N;
corrf(:)=0;
for s=1:(N-1)
    su=0;
    for c=1:(N-s)
        Ic=mean2(data(:,:,c));
        Ics=mean2(data(:,:,c+s));
        for x=1:XS
            for y=1:YS
                su=su+((data(x,y,c)-Ic)*(data(x,y,c+s)-Ics))/Ic/Ics;
            end
        end
    end
    
    corrf(s)=su/(N-s)/XS/YS;
end