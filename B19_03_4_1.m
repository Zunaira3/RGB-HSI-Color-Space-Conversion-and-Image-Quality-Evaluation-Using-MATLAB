%%
%Question 3 implemented on Skyline image
%reading an image
F=imread('skyline.jpg');
F=im2double(F);
%seperating R G B channels
r=F(:,:,1);
g=F(:,:,2);
b=F(:,:,3);
%Converting to HSI using conversion formulas
th=acos((0.5*((r-g)+(r-b)))./((sqrt((r-g).^2+(r-b).*(g-b)))+eps));
H=th;
H(b>g)=2*pi-H(b>g);
H=H/(2*pi);
S=1-3.*(min(min(r,g),b))./(r+g+b+eps);
V=(r+g+b)/3;
HSI=cat(3,H,S,V);
%Applying filter h1 and h2 on the I channel of HSI space
image = im2double(V);
%declaring avg flter 3*3
B = [ 1 1 1;1 1 1; 1 1 1];
%Finding size of image
[r,c] = size (image);
filt_I_h1 = zeros(r,c);
for i = 2:r-1
    for j = 2:c-1
        Im1 = [image(i-1,j-1) image(i-1,j) image(i-1,j+1) ; image(i,j-1) image(i,j) image(i,j+1) ; image(i+1,j-1) image(i+1,j) image(i+1,j+1)];
        IM2 = Im1.*B;
        IM3=IM2/9;
        %avg(i,j)=sum(temp2)/9;
        
        filt_I_h1(i,j) = IM3(1,1) + IM3(1,2) + IM3(1,3) + IM3(2,1) + IM3(2,2) + IM3(2,3) +  IM3(3,1) + IM3(3,2) + IM3(3,3);
        
    end
end
%applying h2 filter on I channel
C = [ 1 4 6 4 1 ;4 16 24 16 4;6 24 36 24 6;4 16 24 16 4; 1 4 6 4 1];
%getting size of image 
[r,c] = size (image);
filt_I_h2 = zeros(r,c);

for i = 3:r-2
    for j = 3:c-2
        V1 = [image(i-2,j-2) image(i-2,j-1) image(i-2,j) image(i-2,j+1) image(i-2,j+2); image(i-1,j-2) image(i-1,j-1) image(i-1,j) image(i-1,j+1) image(i-1,j+2) ; image(i,j-2) image(i,j-1) image(i,j) image(i,j+1) image(i,j+2);image(i+1,j-2) image(i+1,j-1) image(i+1,j) image(i+1,j+1) image(i+1,j+2); image(i+2,j-2) image(i+2,j-1) image(i+2,j) image(i+2,j+1) image(i+2,j+2) ];
        I2 = V1.*C;
        I3=I2/256;

        
        filt_I_h2(i,j) = I3(1,1) + I3(1,2) + I3(1,3)+I3(1,4) + I3(1,5) + I3(2,1) + I3(2,2) + I3(2,3)+I3(2,4) + I3(2,5) +  I3(3,1) + I3(3,2) + I3(3,3)+I3(3,4) + I3(3,5);
    end
end
figure(1)
title('Original image')
subplot 121
imshow(filt_I_h1)
title('h1 filtered  I channel')
subplot 122
imshow(filt_I_h2)
title('h2 filtered   I channel')
%Making an hsi image with H and S same and I as filtered one
hsi=cat(3,H,S,filt_I_h1);
hsi1=cat(3,H,S,filt_I_h2);
%Calling the function which I created for conversion from HSI space to RGB
C=hsitorgb(hsi);
D=hsitorgb(hsi1);

figure(2)
subplot 221
imshow(hsi)
title('HSI Image with I channel Filtered by h1');
subplot 222
imshow(hsi1)
title('HSI Image with I channel Filtered by h2');
subplot 223
imshow(C)
title('RGB Image with I channel Filtered by h1');
subplot 224
imshow(D)
title('RGB Image with I channel Filtered by h2');

[N1,N2] = size (image);
%MSE for h1 filtered image
Im=(F-C).^2;
i1=Im(:,:,1);
i2=Im(:,:,2);
i3=Im(:,:,3);
S1=sum(i1,'all');
S2=sum(i2,'all');
S3=sum(i3,'all');
Avg=(S1+S2+S3)/3;
%S=sum(D(:));
E=1/(N1*N2);
MSE=E.*Avg;

Max=255^2;
disp("PSNR for image with h1 filter :")
PSNR=10*log10(Max/MSE);
disp(PSNR)

% %MSE for h2 filterd image

Im=(F-D).^2;
i1=Im(:,:,1);
i2=Im(:,:,2);
i3=Im(:,:,3);
S1=sum(i1,'all');
S2=sum(i2,'all');
S3=sum(i3,'all');
Avg1=(S1+S2+S3)/3;
E=1/(N1*N2);
MSE1=E.*Avg1;

Max=255^2;
PSNR1=10*log10(Max/MSE1);
disp("PSNR for image with h2 filter :")
disp(PSNR1)

%Function to convert back to RGB space
function C=hsitorgb(hsi)

HV=hsi(:,:,1)*2*pi;

SV=hsi(:,:,2);

IV=hsi(:,:,3);

R=zeros(size(HV));

G=zeros(size(HV));

B=zeros(size(HV));

%RG Sector

id=find((0<=HV)& (HV<2*pi/3));

B(id)=IV(id).*(1-SV(id));

R(id)=IV(id).*(1+SV(id).*cos(HV(id))./cos(pi/3-HV(id)));

G(id)=3*IV(id)-(R(id)+B(id));

%BG Sector

id=find((2*pi/3<=HV)& (HV<4*pi/3));

R(id)=IV(id).*(1-SV(id));

G(id)=IV(id).*(1+SV(id).*cos(HV(id)-2*pi/3)./cos(pi-HV(id)));

B(id)=3*IV(id)-(R(id)+G(id));

%BR Sector

id=find((4*pi/3<=HV)& (HV<2*pi));

G(id)=IV(id).*(1-SV(id));

B(id)=IV(id).*(1+SV(id).*cos(HV(id)-4*pi/3)./cos(5*pi/3-HV(id)));

R(id)=3*IV(id)-(G(id)+B(id));

C=cat(3,R,G,B);

C=max(min(C,1),0);

end














