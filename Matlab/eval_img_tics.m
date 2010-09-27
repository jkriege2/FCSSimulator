% reads images from a diffusion simulation and tries to determine the
% diffusion coefficient using TICS


clear all;


cntmin=20000;   % start from file no <cntmin>
cntmax=21000;   % stop with file no <cntmax>
step=2;        % step through the files by <step> 
sumover=step;   % <sumover> files in the current step will be sumed up to give the next 
                % c(r,dt), so dt = step*dt(Images)

subsize=40;   % size of the region used to calculate the autocorrelation curve in
substep=subsize/2;
               
% Template for the files to read               
filename='./img/test5/img_avg%03d.dat';

% time interval between two subsequent images
dtImages=100e-6; 






dt=dtImages*step;

colors{1}='r';
colors{2}='g';
colors{3}='b';
colors{4}='r:';
colors{5}='g:';
colors{6}='b:';
colors{7}='r--';
colors{8}='g--';
colors{9}='b--';


cntr=0;
rescnt=(cntmax-cntmin)/step;
for cnt=cntmin:step:(cntmax-step)
    cntr=cntr+1;

    % sum over <step> images to generate the next c(t+dt), store c(t) in oldimg
    for s=cnt:(cnt+sumover-1)
        disp([sprintf(filename, s)]);% '   cnt=' num2str(cnt) '   s=' num2str(s)]);
        if (s==cnt)
            img=csvread(sprintf(filename, s));
        else
            img=img+csvread(sprintf(filename, s));
        end
    end
    
    
    imgsize=size(img);
    
    xx=1;
    yy=1;
    for x=1:substep:(imgsize(1)-subsize)
        for y=1:substep:(imgsize(2)-subsize)
            if (cntr==1)
                a=zeros(subsize,subsize,rescnt);
            else
                a=tile{xx,yy};
            end
            %disp(['(x, y)=' num2str([x y])]);
            a(:,:,cntr)=img(x:(x+subsize-1), y:(y+subsize-1));
            tile{xx,yy}=a;
            yy=yy+1;
        end
        xx=xx+1;
        yy=1;
    end
    
end

xx=1;
yy=1;
xcnt=ceil((imgsize(1)-subsize)/substep/2);
ycnt=ceil((imgsize(2)-subsize)/substep/2);
cnt=1;
icnt=1;

for x=1:substep:(imgsize(1)-subsize)
    for y=1:substep:(imgsize(2)-subsize)
        data=tile{xx,yy};
        s=size(data);
        co=tics_corr(data);
        %figure(1)
        if (mod(xx,2)==0 && mod(yy,2)==0) 
            xx/2;
            yy/2;
            xcnt;
            ycnt;
            icnt=(yy/2-1)*xcnt+(xx/2);
            %subplot(xcnt, ycnt, icnt);
            icnt=icnt+1;
        end

        co=co(1:(end-1));
        tau=1:length(co);
        %plot(tau,co,'+');
        
        beta=nlinfit(tau,co,@corrf_diff, [0.1 50 0]);
        
        %plot(tau,co,'+',tau, corrf_diff(beta,tau),'r-');
        c{xx,yy}=co;
        D(xx,yy)=1/4/beta(2);
        cnt=cnt+1;
        figure(2)
        plot(tau,co,'+',tau, corrf_diff(beta,tau),'r-')
        figure(3);
        imagesc(D, [0 20]);
        colorbar;
        yy=yy+1;
    end
    xx=xx+1;
    yy=1;
end

ind=find(abs(D)>=10);
D2=D;
D2(ind)=0;

figure(3);
subplot(2,1,1);
imagesc(D2);
colorbar;
colormap('gray');
subplot(2,1,2);


plot(sum(D2));
