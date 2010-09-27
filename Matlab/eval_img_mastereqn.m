% reads images from a diffusion simulation and tries to determine the
% diffusion coefficient using TICS

if ~exist('eval_several')
    clear all;


    cntmin=1;   % start from file no <cntmin>
    cntmax=20;   % stop with file no <cntmax>
    step=1;         % step through the files by <step> 
    sumover=step;   % <sumover> files in the current step will be sumed up to give the next 
                    % c(r,dt), so dt = step*dt(Images)

    numareas=10; % number of areas with different properties in y direction
    Dareas=[10 5 2.5 1.25 1 0.75 0.5 0.25 0.1 0.05 ]; % diffusion coefficient in the different areas
    
    border=5;

    noiselevel=0;

    % Template for the files to read               
    filename='../img/5/C20_D25pixsteps_01ms_img_cnt%05d';

    % time interval between two subsequent images in seconds
    dtImages=1e-3; 

    % size of one pixel in microns
    pixel_size=200e-3;
end

gaussfilter_onread= exist('gfilter_width') && exist('gfilter_hsize');
if gaussfilter_onread
    if gfilter_width<=0
        gaussfilter_onread=false;
    end
end

temp='';
if exist('filename_post')
    temp=filename_post;
end
filename_post=temp;

dt=sumover*dtImages;

clear img_stack

figure(1)
clf;
figure(2)
clf;
figure(5);
clf;

colors{1}='r';
colors{2}='g';
colors{3}='b';
colors{4}='r:';
colors{5}='g:';
colors{6}='b:';
colors{7}='r--';
colors{8}='g--';
colors{9}='b--';

if exist('makemovie', 'var')
    if makemovie
        try 
            aviobj = avifile([sprintf(filename, cntmax-cntmin+1) '_n' num2str(noiselevel) filename_post '_pic.avi'],'fps',20, 'compression', 'None');%, 'compression', 'Cinepak', 'quality', 100);
        catch MEAVI
            aviobj = close(aviobj);
            aviobj = avifile([sprintf(filename, cntmax-cntmin+1) '_n' num2str(noiselevel)  filename_post '_pic.avi'],'fps',20, 'compression', 'None');%, 'compression', 'Cinepak', 'quality', 100);
        end
    end
end

cntr=0;
rescnt=(cntmax-cntmin)/step;
for cnt=cntmin:step:(cntmax-step)
    cntr=cntr+1;

    % sum over <step> images to generate the next c(t+dt), store c(t) in oldimg
    for s=cnt:(cnt+sumover-1)
        disp([sprintf(filename, s) '.dat']);% '   cnt=' num2str(cnt) '   s=' num2str(s)]);
        if (s==cnt)
            img=csvread([sprintf(filename, s) '.dat']);
        else
            img=img+csvread([sprintf(filename, s) '.dat']);
        end
    end
    
    ma=max(max(img));
    mi=min(min(img));
    
    nimg=img+randn(size(img)).*noiselevel.*abs(ma-mi);
    
    if gaussfilter_onread
        filt=fspecial('gaussian', gfilter_hsize, gfilter_width);
        nimg1=nimg;
        nimg=imfilter(nimg1, filt, 'replicate');
    end
    
    img_stack(:,:,cntr)=nimg;
    figure(1);
    subplot(1,2,1)
    imagesc(img);
    daspect([1 1 1]);
    colorbar;
    subplot(1,2,2);
    imagesc(img_stack(:,:,cntr));
    daspect([1 1 1]);
    colorbar;
    if exist('makemovie', 'var')
        if makemovie
            cmap=colormap('hsv');
            if exist('moviecolormap', 'var')
                cmap=moviecolormap;
            end
            crange=[min(nimg(:)) max(nimg(:))];
            if exist('moviecolorrange', 'var')
                crange=moviecolorrange;
            end
            cimage=(nimg-crange(1))/(crange(2)-crange(1))*(length(cmap)-1)+1;
            err=find(cimage<1.0);
            cimage(err)=1.0;
            err=find(cimage>length(cmap));
            cimage(err)=length(cmap);
            aviobj = addframe(aviobj,im2frame(cimage, cmap));
        end
    end
    
    if cnt<cntmin+12*step
        figure(5);
        set(gcf, 'PaperOrientation', 'landscape');
        subplot(3,4,(cnt-cntmin)/step+1);
        imagesc(img_stack(:,:,cntr));
        daspect([1 1 1]);
        colorbar;
    end
    
end

if exist('makemovie', 'var')
    if makemovie
        aviobj = close(aviobj);
    end
end


h=figure(5);
set(gcf, 'PaperOrientation', 'landscape');
saveas(h, [sprintf(filename, cntmax-cntmin+1) '_n' num2str(noiselevel) filename_post '_pic.fig'])
saveas(h, [sprintf(filename, cntmax-cntmin+1) '_n' num2str(noiselevel) filename_post '_pic.pdf'])

imsize=size(img_stack);
sizex=imsize(1);
sizey=imsize(2);
sizet=imsize(3);


clear p_res
clear p_resf
clear residuum

gauss5 = [0.0005    0.0050    0.0109    0.0050    0.0005;    0.0050    0.0522    0.1141    0.0522    0.0050;    0.0109    0.1141    0.2491    0.1141    0.0109;    0.0050    0.0522    0.1141    0.0522    0.0050;    0.0005    0.0050    0.0109    0.0050    0.0005];
tic
for x=border:(sizex-border)
    for y=border:(sizey-border)
        dataslice=img_stack((x-1):(x+1),(y-1):(y+1),:);
        %f = @(p)mastereqn_optfun(p, dataslice);
        flsq = @(p)mastereqn_optfun_lsq(p, dataslice);
        %flsq = @(p)mastereqn_optfun1_lsq(p, dataslice, dt);
        p0=1;
        options = optimset('LargeScale','off');
        %[p,fval,exitflag,output] = fminunc(f,p0,options);
        [p,r]=lsqnonlin(flsq,p0);
        p_res(x-border+1,y-border+1)=p;%*pixel_size^2;
        residuum(x-border+1,y-border+1)=r;
        try
            figure(2);
            subplot(2,3,1);
            imagesc(p_res);
            daspect([1 1 1]);
            %colorbar
            subplot(2,3,2);
            imagesc(residuum);
            daspect([1 1 1]);
            %colorbar
            subplot(2,3,4);
            mp_res=mean(p_res,1);
            plot(mp_res);
            xlim([1 sizey-2*border]);
            ylim([min([0 min(mp_res)]) 1.1*max(mp_res(:))]);
            meanv=0;
            stddev=0;
            if exist('numareas')
                areawidthbase=sizey/numareas;
                areawidthshort=areawidthbase-border;
                yy=1;
                for yc=1:numareas%areawidth:(sizey-2*border)
                    areawidth=areawidthbase;
                    if yc==1
                        areawidth=areawidthshort;
                    end
                    sizpr=size(p_res);
                    sizpy=sizpr(2);
                    meana=p_res(:,min(yy,sizpy):min(sizpy,(yy+areawidth)));
                    meanv(yc)=mean(meana(:));
                    stddev(yc)=std(meana(:));
                    line([yy (yy+areawidth)], [meanv(yc) meanv(yc)], 'Color', 'r');
                    line([yy (yy+areawidth)], [meanv(yc) meanv(yc)]-stddev(yc), 'Color', 'm');
                    line([yy (yy+areawidth)], [meanv(yc) meanv(yc)]+stddev(yc), 'Color', 'm');
                    yy=yy+areawidth;
                end
                if mod(y, areawidthbase)==0
                    subplot(2,3,5);
                    p0=[0 1];
                    p = lsqcurvefit(@fit_linear,p0,Dareas,meanv);
                    p1 = lsqcurvefit(@fit_linear,p0,Dareas(1:4),meanv(1:4));
                    p2 = lsqcurvefit(@fit_linear,p0,Dareas(end-4:end),meanv(end-4:end));
                    errorbar(Dareas, meanv, stddev, 'xb');
                    hold on
                    xlabel('supplied values');
                    ylabel('estimated');
                    limits=ylim;
                    ylim([0 limits(2)]);
                    limits=xlim;
                    xlim([0 limits(2)]);
                    xl=0:(limits(2)/5):limits(2);
                    plot(xl, fit_linear(p, xl), 'r-');
                    plot(xl, fit_linear(p1, xl), 'g--');
                    plot(xl, fit_linear(p2, xl), 'b--');
                    text(0.1,0.85, ['f(x) = ' num2str(p(1)) ' + ' num2str(p(2)) '* x'], 'Units', 'Normalized');
                    text(0.1,0.75, ['f(x) = ' num2str(p1(1)) ' + ' num2str(p1(2)) '* x'], 'Units', 'Normalized');
                    text(0.1,0.65, ['f(x) = ' num2str(p2(1)) ' + ' num2str(p2(2)) '* x'], 'Units', 'Normalized');
                    hold off;
                end
            else
                meanv=mean(p_res(:));
                stddev=std(p_res(:));
                line([1 sizey-2*border], [meanv meanv], 'Color', 'r');
                line([1 sizey-2*border], [meanv meanv]-stddev, 'Color', 'm');
                line([1 sizey-2*border], [meanv meanv]+stddev, 'Color', 'm');
                text(0.5,0.8, ['mean=' num2str(meanv)], 'Units', 'Normalized');
                text(0.5,0.65, ['stddev=' num2str(stddev)], 'Units', 'Normalized');
            end
            p_resf=imfilter(p_res, fspecial('gaussian', [5 5], 0.8));
            %p_resf=conv2(p_res, gauss5);
            subplot(2,3,3);
            imagesc(p_resf);
            daspect([1 1 1]);
            subplot(2,3,6);
            %plot(mean(p_resf(5:(end-5),5:(end-5)),1));
            plot(mean(p_resf,1));
            %xlim([1 sizex-2*border]);
            xlim([1 sizey-10-2*border]);
        catch ME
            disp('ERROR:');
            disp(ME.message);
        end
    end
end
toc
meanv=mean(p_res(:))
stddev=std(p_res(:))

try
    figure(2);
    set(gcf, 'PaperOrientation', 'landscape');
    subplot(2,3,1);
    imagesc(p_res);
    daspect([1 1 1]);
    colorbar;
    subplot(2,3,4);
    mp_res=mean(p_res,1);
    plot(mp_res);
    xlim([1 sizey-2*border]);
    ylim([min([0 min(mp_res)]) 1.1*max(mp_res)]);
%     line([1 sizey-2*border], [meanv meanv]);
%     line([1 sizey-2*border], [meanv meanv]-stddev);
%     line([1 sizey-2*border], [meanv meanv]+stddev);
%     text(0.5,0.8, ['mean=' num2str(meanv)], 'Units', 'Normalized');
%     text(0.5,0.65, ['stddev=' num2str(stddev)], 'Units', 'Normalized');
    if exist('numareas')
        areawidthbase=sizey/numareas;
        areawidthshort=areawidthbase-border;
        yy=1;
        for yc=1:numareas%areawidth:(sizey-2*border)
            areawidth=areawidthbase;
            if yc==1
                areawidth=areawidthshort;
            end
            sizpr=size(p_res);
            sizpy=sizpr(2);
            meana=p_res(:,min(yy,sizpy):min(sizpy,(yy+areawidth)));
            meanv(yc)=mean(meana(:));
            stddev(yc)=std(meana(:));
            line([yy (yy+areawidth)], [meanv(yc) meanv(yc)], 'Color', 'r');
            line([yy (yy+areawidth)], [meanv(yc) meanv(yc)]-stddev(yc), 'Color', 'm');
            line([yy (yy+areawidth)], [meanv(yc) meanv(yc)]+stddev(yc), 'Color', 'm');
            yy=yy+areawidth;
        end
        subplot(2,3,5);
        p0=[0 1];
        p = lsqcurvefit(@fit_linear,p0,Dareas,meanv);
        plot(Dareas, meanv)
        errorbar(Dareas, meanv, stddev, 'xb');
        hold on
        xlabel('supplied values');
        ylabel('estimated');
        limits=ylim;
        ylim([0 limits(2)]);
        limits=xlim;
        xlim([0 limits(2)]);
        xl=0:(limits(2)/5):limits(2);
        plot(xl, fit_linear(p, xl), 'r-');
        text(0.1,0.85, ['f(x) = ' num2str(p(1)) ' + ' num2str(p(2)) '* x'], 'Units', 'Normalized');
        hold off;
        A=transpose([Dareas; meanv; stddev]);
        dlmwrite([sprintf(filename, cntmax-cntmin+1) '_n' num2str(noiselevel) filename_post '_result.txt'], A, ';');
    else
        meanv=mean(p_res(:));
        stddev=std(p_res(:));
        line([1 sizey-2*border], [meanv meanv], 'Color', 'r');
        line([1 sizey-2*border], [meanv meanv]-stddev, 'Color', 'm');
        line([1 sizey-2*border], [meanv meanv]+stddev, 'Color', 'm');
        text(0.5,0.8, ['mean=' num2str(meanv)], 'Units', 'Normalized');
        text(0.5,0.65, ['stddev=' num2str(stddev)], 'Units', 'Normalized');
    end


    subplot(2,3,2);
    imagesc(residuum);
    daspect([1 1 1]);
    colorbar;

     p_resf=imfilter(p_res, fspecial('gaussian', [5 5], 0.8));
      %p_resf=conv2(p_res, gauss5);

     subplot(2,3,3);
     imagesc(p_resf);
     daspect([1 1 1]);
     colorbar;
     subplot(2,3,6);
     %plot(mean(p_resf(5:(end-5),5:(end-5)),1));
     plot(mean(p_resf,1));
     %xlim([1 sizex-10-2*border]);
     xlim([1 sizey-2*border]);
     ylim([min([0 min(mean(p_resf,1))]) 1.1*max(mean(p_resf,1))]);


    saveas(gcf, [sprintf(filename, cntmax-cntmin+1) '_n' num2str(noiselevel) filename_post '.fig'])
    saveas(gcf, [sprintf(filename, cntmax-cntmin+1) '_n' num2str(noiselevel) filename_post '.pdf'])
catch ME
    disp('ERROR:');
    disp(ME.message);
end
