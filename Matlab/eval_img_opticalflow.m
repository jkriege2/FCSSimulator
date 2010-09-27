% reads images from a diffusion simulation and tries to determine the
% diffusion coefficient using TICS


clear all;


cntmin=10000;   % start from file no <cntmin>
cntmax=10050;   % stop with file no <cntmax>
step=1;        % step through the files by <step> 
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
    
    img_stack(:,:,cntr)=img;
    
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

figure(1)
implay(img_stack);

[Vx Vy]=OpticalFlow(img_stack, 5,10);

figure(2);
subplot(4,1,1);
imagesc(Vx);
subplot(4,1,2);
imagesc(Vy);
subplot(4,1,3);
Vabs=sqrt(Vx.^2+Vy.^2);
imagesc(Vabs);
subplot(4,1,4);
meanfilter=[1 2 4 2 1];
meanfilter=meanfilter./sum(meanfilter);
plot(filter(meanfilter,1,mean(Vabs,1)));
hold on
plot(mean(Vabs,1),'r');
hold off;


