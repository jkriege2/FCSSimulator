clear all;
cntmin=1000 ;
cntmax=20000;
step=1000;
filename='./img/test5/img_avg%03d.dat';
% border=15;
% segwidth=180/4;
% sumover=1;
% D=[10 100 1000 10000];
% dt=100e-6;
% ranges=4;
border=10;
segwidth=200/5;
sumover=10;
D=[1000 5000 10000 15000 20000];
dt=100e-6;
ranges=5;


step=max(step,sumover);
colors{1}='r';
colors{2}='g';
colors{3}='b';
colors{4}='r:';
colors{5}='g:';
colors{6}='b:';
colors{7}='r--';
colors{8}='g--';
colors{9}='b--';

for c=1:ranges
    legendtext{c}=[num2str(D(c)) ' {\mu}m^2/s'];
end

cntr=0;
for cnt=cntmin:step:(cntmax-step)
    cntr=cntr+1;
    for s=cnt:(cnt+sumover-1)
        disp([sprintf(filename, s)]);% '   cnt=' num2str(cnt) '   s=' num2str(s)]);
        if (s==cnt)
            img=csvread(sprintf(filename, s));
        else
            img=img+csvread(sprintf(filename, s));
        end
    end
    figure(1);
    subplot(2,1,1)
    imagesc(img);
    colorbar;
    subplot(2,1,2)
    plot(1:length(img), sum(img));
    pause(0.01)
    %if (cnt==0) 
    %    img1=img;
    %else
    %    img1=img1+img;
    %end
    %figure(2)
    %imagesc(img1);
    
    img1=zeros(length(img),ranges*(segwidth-2*border));
    
    for c=1:ranges
        %disp([ '[' num2str(((c-1)*segwidth+1+border)) ' .. ' num2str((c*segwidth-border)) ']']);
        data=img(:,((c-1)*segwidth+1+border):(c*segwidth-border));
        a=acf(data);
        result(cntr,c,1)=mean2(data);
        result(cntr,c,2)=std2(data);
        img1(:,((c-1)*(segwidth-2*border)+1):(c*(segwidth-2*border)))=a;
        cut=a(length(a)/2+1,:);
        beta=nlinfit(1:length(cut),cut,@gauss_bias, [mean(cut) max(cut) length(cut)/2 1]);
        result(cntr,c,3:6)=beta;
        
        if cnt==cntmin
            avgacf{c}=a;
        else 
            avgacf{c}=avgacf{c}+a;
        end
        avgimg1(:,((c-1)*(segwidth-2*border)+1):(c*(segwidth-2*border)))=avgacf{c};

        figure(4)
        plot(1:cntr, result(1:cntr,c,1),colors{c});
        hold on;
    end
    hold off;
    
    figure(2)
    subplot(2,2,1);
    imagesc(img1);
    colorbar;
    subplot(2,2,2);
    imagesc(avgimg1./cntr);
    colorbar;
    subplot(2,2,3);
    s=size(img1);
    completelen=s(2);
    singlelen=completelen/ranges;
    plot(1:completelen, img1(length(img1)/2+1,:), 'b+');
    hold on;
    for c=1:ranges
        plot(((c-1)*singlelen+1):0.1:(c*singlelen), gauss_bias(result(cntr,c,3:6), 1:0.1:singlelen), 'r');
    end
    hold off;
    subplot(2,2,4);
    for c=1:ranges
        s=size(avgimg1);
        completelen=s(2);
        singlelen=completelen/ranges;
        a=avgacf{c}/cntr;
        cut=a(length(a)/2+1,:);
        plot(((c-1)*singlelen+1):1:(c*singlelen),cut, 'b+');
        hold on;
    end
    hold off;
    
    
end

w=[];
for c=1:ranges
    figure(3)
    subplot(2,1,1)
    plot(1:cntr, abs(result(1:cntr,c,6)),colors{c});
    w(c)=mean(abs(result(1:cntr,c,6)));
    hold on;

    a=avgacf{c}/cntr;
    cut=a(length(a)/2+1,:);
    beta=nlinfit(1:length(cut),cut,@gauss_bias, [mean(cut) max(cut) length(cut)/2 1]);
    aw(c)=abs(beta(4));

    
    figure(2)
    subplot(2,2,2);
    imagesc(avgimg1./cntr);
    colorbar;
    subplot(2,2,4);
    s=size(avgimg1);
    completelen=s(2);
    singlelen=completelen/ranges;
    plot(((c-1)*singlelen+1):1:(c*singlelen),cut, 'b+');
    hold on;
    plot(((c-1)*singlelen+1):0.1:(c*singlelen), gauss_bias(beta, 1:0.1:singlelen), 'r');
    
end
figure(2)
hold off;
figure(3);
hold off;
legend(legendtext(1:ranges));
title('ACF width');
subplot(2,1,2);
p=polyfit(sqrt(D),w,1);
Dvec=D(1):((D(end)-D(1))/100):D(end);
disp(['fit/sqrt(6t) = ' num2str(p(1)/sqrt(6*dt*sumover)) ]);
plot(sqrt(D),w,'+',sqrt(D),aw,'o',sqrt(Dvec),polyval(p,sqrt(Dvec)),'-');%, sqrt(Dvec), sqrt(Dvec*6*dt*sumover), ':');
title('ACF width vs. sqrt(D)');
xlabel('sqrt(D)');
ylabel('ACF width [pixel]');



figure(4)
for c=1:ranges
    plot(1:cntr, result(1:cntr,c,1),colors{c});
    hold on;
end
hold off;
legend(legendtext(1:ranges));
title('average photon count');
%figure(2)
%imagesc(img1);