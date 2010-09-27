% reads images from a diffusion simulation and tries to determine the
% diffusion coefficient using TICS


clear all;


cntmin=20000;   % start from file no <cntmin>
cntmax=20051;   % stop with file no <cntmax>
step=1;         % step through the files by <step> 
sumover=step;   % <sumover> files in the current step will be sumed up to give the next 
                % c(r,dt), so dt = step*dt(Images)

subsize=44;   % size of the region used to calculate the autocorrelation curve in
substep=subsize/2;

border=5;
               
% Template for the files to read               
filename='./img/test5/img_avg%03d.dat';

% time interval between two subsequent images
dtImages=1; 

dt=sumover*dtImages;





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
    
end

imsize=size(img_stack);
sizex=imsize(1);
sizey=imsize(2);
sizet=imsize(3);

cx=0;
for x=border:(sizex-border)
    cx=cx+1;
end

cy=0;
for y=border:(sizey-border)
    cy=cy+1;
end

param_count=cx*cy;

tic

flsq = @(p)mastereqn_optfun_complete_lsq(p, img_stack, dt, border);
p0(1:param_count)=1;
options = optimset('LargeScale','on');
p=lsqnonlin(flsq,p0);
p_res=reshape(p, cx,cy);
toc

figure(2);
subplot(2,2,1);
imagesc(p_res);
subplot(2,2,3);
plot(mean(p_res,1));

p_resf=imfilter(p_res, fspecial('gaussian', [5 5], 0.8));

subplot(2,2,2);
imagesc(p_resf);
subplot(2,2,4);
plot(mean(p_resf,1));

