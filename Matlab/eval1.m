figure(1)
clear all;
icnt=1;
icnt1=6;
for cnt=1:10
    fn=sprintf('10ms/img%03d.dat', cnt-1);
    img{cnt}=csvread(fn);
    im=img{cnt};
    figure(1);
    subplot(4,5,icnt);
    imagesc(im);%, [0 30000]);
    axis off;
    a=im;
    a(:,1:end/2)=acf(im(:,1:end/2));
    a(:,1:end/2)=a(:,1:end/2)-mean(mean(a(1:5,1:5)));
    disp(mean(mean(a(1:5,1:5))));
    a(:,1:end/2)=a(:,1:end/2)./max(max(a(:,1:end/2)));
    a(:,end/2+1:end)=acf(im(:,end/2+1:end));
    a(:,end/2+1:end)=a(:,end/2+1:end)-mean(mean(a(1:5, end-5:end)));
    disp(mean(mean(a(1:5, end-5:end))));
    a(:,end/2+1:end)=a(:,end/2+1:end)./max(max(a(:,end/2+1:end)));
    figure(1);
    subplot(4,5,icnt1);
    imagesc(a);
    corr{cnt}=a;
    axis off;
    icnt=icnt+1;
    icnt1=icnt1+1;
    if (icnt==6) 
        icnt=11;
        icnt1=16;
    end
end

figure(2)
c=corr{1};
for cnt=2:10
    c=c+corr{cnt};
end
c=c./10;
imagesc(c);
figure(3);
c1=c(end/2,1:end/2);
c2=c(end/2,(end/2+1):end);
b1=nlinfit(1:length(c1), c1, @gauss, [length(c1)/2 1 1]);
b1(2)*0.1
b2=nlinfit(1:length(c2), c2, @gauss, [length(c1)/2 1 1]);
b2(2)*0.1
plot(c1, 'r+');
hold on
plot(1:0.1:length(c1), gauss(b1, 1:0.1:length(c1)), 'r');
plot(c2, 'bo');
plot(1:0.1:length(c2), gauss(b2, 1:0.1:length(c2)), 'b');
hold off
figure(4)
subplot(1,2,1)
surf(c(1:end,1:end/2))
subplot(1,2,2)
surf(c(1:end,end/2:end))
% for cnt=1:2
%     fn=sprintf('traj%03d.dat', cnt-1);
%     traj{cnt}=csvread(fn);
%     t=traj{cnt};
%     figure(2);
%     subplot(2,1,cnt);
%     plot3(t(:,1), t(:,2), t(:,3));
%     axis on;
% end