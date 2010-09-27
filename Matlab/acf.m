function a = acf(img)
s=size(img);
mask=img;
sigma=0.3;
for xx=1:s(1)
    for yy=1:s(2)
        mask(xx,yy)=sin(xx./s(1).*pi)*sin(yy./s(2).*pi);
        %mask(xx,yy)=0.5*(1-cos(2*pi*xx/(s(1)-1)))*0.5*(1-cos(2*pi*yy/(s(2)-1)));
        %mask(xx,yy)=exp(-0.5*(xx-(s(1)-1)/2)^2/(sigma*(s(1)-1)/2)^2-0.5*(yy-(s(2)-1)/2)^2/(sigma*(s(2)-1)/2)^2);
    end
end
%figure(10);
%imagesc(mask);
img1=img;
%img1=img1.*mask;
meanI=mean(mean(img1));
R1=fft2(img1);
a=fftshift(real(ifft2(abs(R1).^2))./(meanI.^2))-1;

%a=conv2(img, img)./(meanI.^2)-1;