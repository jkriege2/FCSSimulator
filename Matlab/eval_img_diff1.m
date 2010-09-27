% reads images from a diffusion simulation and tries to determine the
% diffusion coefficient at every position by calculating:
%
%        laplace(c)
%  D = --------------
%         dc / dt
%
% which may directly be obtained from the diffusion equation. THis is
% similar to what is done with the law of intensity conservation for image
% flow estimation
%
% The concentrations c(r,t) are obtained by averaging over different
% numbers of "noise images".
%
% The laplacian is estimated by a laplace filter and the time derivative by
% 
%            c(t+dt) - c(t)
%  dc/dt = ------------------
%                  dt


clear all;


cntmin=20000;   % start from file no <cntmin>
cntmax=21000;  % stop with file no <cntmax>
step=10;     % step through the files by <step> ... note that all <step> 
               % files in the current step will be sumed up to give the next 
               % c(r,dt), so dt = step*dt(Images)
               
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


img=csvread(sprintf(filename, cntmin));
img(:,:)=0;
oldimg=img;

cntr=0;
rescnt=(cntmax-cntmin)/step;
for cnt=cntmin:step:(cntmax-step)
    cntr=cntr+1;

    % sum over <step> images to generate the next c(t+dt), store c(t) in oldimg
    oldoldimg=oldimg;
    oldimg=img;
    for s=cnt:(cnt+step-1)
        disp([sprintf(filename, s)]);% '   cnt=' num2str(cnt) '   s=' num2str(s)]);
        img=img+csvread(sprintf(filename, s));
    end
    
    if (cntr>=3)
    
        aoldimg=oldimg/((cntr-1)*step);
        aoldoldimg=oldoldimg/((cntr-2)*step);
        aimg=img/(cntr*step);
        figure(1);
        subplot(2,3,1)
        imagesc(aoldoldimg);
        colorbar;
        title('c(x,y,t)');
        subplot(2,3,4)
        plot(1:length(aoldimg), aoldimg(end/2,:));
        title('c(x,y,t)');
        
        subplot(2,3,2)
        imagesc(aoldimg);
        colorbar;
        title('c(x,y,t+dt)');
        subplot(2,3,5)
        plot(1:length(aoldimg), aoldimg(end/2,:));
        title('c(x,y,t+dt)');
        
        subplot(2,3,3)
        imagesc(aimg);
        colorbar;
        title('c(x,y,t+2*dt)');
        subplot(2,3,6)
        plot(1:length(aimg), aimg(end/2,:));
        title('c(x,y,t+2*dt)');
        
        dcdt=(imfilter(aimg,fspecial('gaussian'))-imfilter(aoldoldimg,fspecial('gaussian')))./(2*dt);
        laplaceC=imfilter(aimg, fspecial('laplacian'));%del2(aimg);
        
        dcdt=imresize(dcdt, 0.5);
        laplaceC=imresize(laplaceC, 0.5);
        
        figure(3);
        subplot(2,1,1);
        imagesc(dcdt);
        title('dc/dt');
        colorbar;
        subplot(2,1,2);
        imagesc(laplaceC);
        colorbar;
        title('Laplace(c(t+dt))');
        
        D=imfilter(dcdt./laplaceC, fspecial('gaussian'));
        
        if (cntr==3)
            Da=D;
        else
            Da=Da+D;
        end
        mean(mean(D))
        std2(D)
        figure(2);
        subplot(ceil(sqrt(rescnt)), ceil(sqrt(rescnt)), cntr);
        imagesc(D, [0 100]);
        colorbar;
        figure(4);
        subplot(ceil(sqrt(rescnt)), ceil(sqrt(rescnt)), cntr);
        plot(1:length(D), mean(D(:,:)));
        
        figure(5)
        subplot(2,1,1);
        subplot(2,1,1);
        imagesc(Da/(cntr-2),[0 100]);
        colorbar
        subplot(2,1,2);
        DD=Da/(cntr-2);
        DD(find(DD<0))=0;
        DD(find(DD>100))=100;
        plot(sum(DD));
    end
end
figure(5)
subplot(2,1,1);
imagesc(Da/(cntr-2),[0 100]);
colorbar
subplot(2,1,2);
DD=Da/(cntr-2);
DD(find(DD<0))=0;
DD(find(DD>100))=100;
plot(sum(DD));