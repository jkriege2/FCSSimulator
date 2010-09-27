% taken from:
% http://chi3x10.wordpress.com/2008/06/09/optical-flow-in-matlab/

function [Vx,Vy] = OpticalFlow(images,alpha,iterations)  
 %// Calculating optical flow of a sequence of images.  
 %// images : 3D array that contains a sequence of images. size of images is (imageHeight, imageWidth, frame number)  
 %// alpha  
 %// iterations.  
 [height,width,frames]=size(images);  
   
 %//initialzation of u and v  
 Vx = zeros(height,width);  
 Vy = zeros(height,width);  
   
 for k = 1:frames-1  
   
     % //initialization of Ex Ey and Et  
     Ex = zeros(height-1,width-1,frames-1);  
     Ey = zeros(height-1,width-1,frames-1);  
     Et = zeros(height-1,width-1,frames-1);  
   
     %//calculating Ex Ey and Et in frame k.  
     for x = 2:width-1  
         for y = 2:height-1  
             Ex(y,x,k) = (images(y+1,x+1,k)-images(y+1,x,k)+images(y,x+1,k)...  
                 -images(y,x,k)+images(y+1,x+1,k+1)-images(y+1,x,k+1)...  
                 +images(y,x+1,k+1)-images(y,x,k+1))/4;  
   
             Ey(y,x,k) = (images(y,x,k)-images(y+1,x,k)+images(y,x+1,k)...  
                 -images(y+1,x+1,k)+images(y,x,k+1)-images(y+1,x,k+1)...  
                 +images(y,x+1,k+1)-images(y+1,x+1,k+1))/4;  
   
             Et(y,x,k) = (images(y+1,x,k+1)-images(y+1,x,k)+images(y,x,k+1)...  
                 -images(y,x,k)+images(y+1,x+1,k+1)-images(y+1,x+1,k)...  
                 +images(y,x+1,k+1)-images(y,x+1,k))/4;  
         end  
     end  
   
     for nn = 1:iterations  
         for x = 2:width-1  
             for y = 2:height-1  
   
                 Vxbar = (Vx(y-1,x)+Vx(y,x+1)+Vx(y+1,x)+Vx(y,x-1))/6+...  
                      (Vx(y-1,x-1)+Vx(y-1,x+1)+Vx(y+1,x+1)+Vx(y+1,x-1))/12;  
   
                 Vybar = (Vy(y-1,x)+Vy(y,x+1)+Vy(y+1,x)+Vy(y,x-1))/6+...  
                     (Vy(y-1,x-1)+Vy(y-1,x+1)+Vy(y+1,x+1)+Vy(y+1,x-1))/12;  
   
                 %// chapter 12 of Horn's paper  
                 temp = (Ex(y,x,k)*Vxbar+Ey(y,x,k)*Vybar+Et(y,x,k))/(alpha^2 + Ex(y,x,k)^2 + Ey(y,x,k)^2);  
                 %// update u and v  
                 Vx(y,x) = Vxbar-Ex(y,x,k)*temp;  
                 Vy(y,x) = Vybar-Ey(y,x,k)*temp;  
             end  
         end  
     end  
   
 end 