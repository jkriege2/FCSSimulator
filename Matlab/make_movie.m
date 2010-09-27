% reads images from a diffusion simulation and tries to determine the
% diffusion coefficient using TICS


clear all;


cntmin=0;   % start from file no <cntmin>
cntmax=233;   % stop with file no <cntmax>
step=1;         % step through the files by <step> 
sumover=step;   % <sumover> files in the current step will be sumed up to give the next 
                % c(r,dt), so dt = step*dt(Images)


% Template for the files to read               
filename='../img/3/D10_01ms_img_avg%05d.dat';
aviobj = avifile('../3_01ms.avi','fps',15, 'compression', 'None');%, 'compression', 'Cinepak', 'quality', 100);
cm=jet;

cntr=0;

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
    imagesc(img);
    colorbar
    aviobj = addframe(aviobj,gca);
    
end
%close(fig);
%movie(img_stack)
imsize=size(img_stack);
sizex=imsize(1);
sizey=imsize(2);
sizet=imsize(3);

aviobj = close(aviobj);

