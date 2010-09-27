clear all;

cntmin=1;   % start from file no <cntmin>
cntmax=100;   % stop with file no <cntmax>
step=1;         % step through the files by <step> 
sumover=step;   % <sumover> files in the current step will be sumed up to give the next 
                % c(r,dt), so dt = step*dt(Images)
border=5;
noiselevel=0;
% size of one pixel in microns
pixel_size=200e-3;

numareas=10; % number of areas with different properties in y direction
Dareas=[ 10 5 2.5 1.25 1 0.75 0.5 0.25 0.1 0.05  ]; % diffusion coefficient in the different areas
countmax=[20 100 1000];

dtImages=1e-3;

moviecolormap=colormap('hsv');
moviecolorrange=[1e4 16e4];

filename_post='_gaussianfiltered_w1';

gfilter_width=1;
gfilter_hsize=[3 3];

eval_several=true;

for nl=1:2
    noiselevel=nl*0.02;
    for cnt=1:length(countmax)
        filename='../img/5/C20_D25pixsteps_01ms_img_cnt%05d';
        cntmax=countmax(cnt);
        makemovie=true;
        run eval_img_mastereqn 

    end
end

