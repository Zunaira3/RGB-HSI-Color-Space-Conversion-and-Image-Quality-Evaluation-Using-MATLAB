
clear all
close all
%Reading Image 
A=imread('cameraman.png');
A=im2double(A);
%Seperating red green  and blue channels
Red=A(:,:,1);
Green=A(:,:,2);
Blue=A(:,:,3);
%plotting image 
figure (1)
subplot 221
imshow(A)
title('Cameraman Image')
%Displaying red channel
subplot 222
imshow(Red)
title('Red Channel')
%Displaying green channel
subplot 223
imshow(Green)
title('Green Channel')
%Displaying blue channel
subplot 224
imshow(Blue)
title('Blue Channel')

image = im2double(Red);
%declaring avg flter 3*3
B = [ 1 1 1;1 1 1; 1 1 1];
%Finding size of image
[r,c] = size (image);
IMG = zeros(r,c);
for i = 2:r-1
    for j = 2:c-1
        Im1 = [image(i-1,j-1) image(i-1,j) image(i-1,j+1) ; image(i,j-1) image(i,j) image(i,j+1) ; image(i+1,j-1) image(i+1,j) image(i+1,j+1)];
        IM2 = Im1.*B;
        IM3=IM2/9;
        IMG(i,j) = IM3(1,1) + IM3(1,2) + IM3(1,3) + IM3(2,1) + IM3(2,2) + IM3(2,3) +  IM3(3,1) + IM3(3,2) + IM3(3,3);
        
    end
end
image1 = im2double(Green);
%declaring avg flter 3*3
B = [ 1 1 1;1 1 1; 1 1 1];
%Finding size of image
[r,c] = size (image1);
IMG1 = zeros(r,c);
for i = 2:r-1
    for j = 2:c-1
        Im1 = [image1(i-1,j-1) image1(i-1,j) image1(i-1,j+1) ; image1(i,j-1) image1(i,j) image1(i,j+1) ; image1(i+1,j-1) image1(i+1,j) image1(i+1,j+1)];
        IM2 = Im1.*B;
        IM3=IM2/9;
        IMG1(i,j) = IM3(1,1) + IM3(1,2) + IM3(1,3) + IM3(2,1) + IM3(2,2) + IM3(2,3) +  IM3(3,1) + IM3(3,2) + IM3(3,3);
        
    end
end

image2 = im2double(Blue);
%declaring avg flter 3*3
B = [ 1 1 1;1 1 1; 1 1 1];
%Finding size of image
[r,c] = size (image2);
IMG2 = zeros(r,c);
for i = 2:r-1
    for j = 2:c-1
        Im1 = [image2(i-1,j-1) image2(i-1,j) image2(i-1,j+1) ; image2(i,j-1) image2(i,j) image2(i,j+1) ; image2(i+1,j-1) image2(i+1,j) image2(i+1,j+1)];
        IM2 = Im1.*B;
        IM3=IM2/9;
        IMG2(i,j) = IM3(1,1) + IM3(1,2) + IM3(1,3) + IM3(2,1) + IM3(2,2) + IM3(2,3) +  IM3(3,1) + IM3(3,2) + IM3(3,3);
        
    end
end

h1=cat(3,IMG,IMG1,IMG2);
C = [ 1 4 6 4 1 ;4 16 24 16 4;6 24 36 24 6;4 16 24 16 4; 1 4 6 4 1];
[r,c] = size (Red);
Img = zeros(r,c);
Im = im2double(Red);
%applying filter on image without using function.It will apply filter one
%by one on each pixel and give output
for i = 3:r-2
    for j = 3:c-2
        i1 = [Im(i-2,j-2) Im(i-2,j-1) Im(i-2,j) Im(i-2,j+1) Im(i-2,j+2); Im(i-1,j-2) Im(i-1,j-1) Im(i-1,j) Im(i-1,j+1) Im(i-1,j+2) ; Im(i,j-2) Im(i,j-1) Im(i,j) Im(i,j+1) Im(i,j+2);Im(i+1,j-2) Im(i+1,j-1) Im(i+1,j) Im(i+1,j+1) Im(i+1,j+2); Im(i+2,j-2) Im(i+2,j-1) Im(i+2,j) Im(i+2,j+1) Im(i+2,j+2) ];
        i2 = i1.*C;
        i3=i2/256;
        Img(i,j) = i3(1,1) + i3(1,2) + i3(1,3)+i3(1,4) + i3(1,5) + i3(2,1) + i3(2,2) + i3(2,3)+i3(2,4) + i3(2,5) +  i3(3,1) + i3(3,2) + i3(3,3)+i3(3,4) + i3(3,5);
    end
end

[m,n] = size (Green);
X = zeros(m,n);
X1 = im2double(Green);
for i = 3:r-2
    for j = 3:c-2
        Y = [X1(i-2,j-2) X1(i-2,j-1) X1(i-2,j) X1(i-2,j+1) X1(i-2,j+2); X1(i-1,j-2) X1(i-1,j-1) X1(i-1,j) X1(i-1,j+1) X1(i-1,j+2) ; X1(i,j-2) X1(i,j-1) X1(i,j) X1(i,j+1) X1(i,j+2);X1(i+1,j-2) X1(i+1,j-1) X1(i+1,j) X1(i+1,j+1) X1(i+1,j+2); X1(i+2,j-2) X1(i+2,j-1) X1(i+2,j) X1(i+2,j+1) X1(i+2,j+2) ];
        X2 = Y.*C;
        X3=X2/256;
        X(i,j) = X3(1,1) + X3(1,2) + X3(1,3)+X3(1,4) + X3(1,5) + X3(2,1) + X3(2,2) + X3(2,3)+X3(2,4) + X3(2,5) +  X3(3,1) + X3(3,2) + X3(3,3)+X3(3,4) + X3(3,5);
    end
end


[y,z] = size (Blue);
S = zeros(y,z);
S1 = im2double(Blue);
for i = 3:r-2
    for j = 3:c-2
        Y1 = [S1(i-2,j-2) S1(i-2,j-1) S1(i-2,j) S1(i-2,j+1) S1(i-2,j+2); S1(i-1,j-2) S1(i-1,j-1) S1(i-1,j) S1(i-1,j+1) S1(i-1,j+2) ; S1(i,j-2) S1(i,j-1) S1(i,j) S1(i,j+1) S1(i,j+2);S1(i+1,j-2) S1(i+1,j-1) S1(i+1,j) S1(i+1,j+1) S1(i+1,j+2); S1(i+2,j-2) S1(i+2,j-1) S1(i+2,j) S1(i+2,j+1) S1(i+2,j+2) ];
        S2 = Y1.*C;
        S3=S2/256;
        S(i,j) = S3(1,1) + S3(1,2) + S3(1,3)+S3(1,4) + S3(1,5) + S3(2,1) + S3(2,2) + S3(2,3)+S3(2,4) + S3(2,5) +  S3(3,1) + S3(3,2) + S3(3,3)+S3(3,4) + S3(3,5);
    end
end
h2=cat(3,Img,X,S);
figure(2)
subplot 131
imshow(A)
title('Original image')
subplot 132
imshow(h1)
title('After applying h1 filter')
subplot 133
imshow(h2)
title('After applying h2 filter')

[N1, N2]=size(A);
%MSE for h1 filtered image
Im=(A-h1).^2;
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
disp("MSE of image with h1 filter :")
disp(MSE)
Max=255^2;
disp("PSNR for image with h1 filter :")
PSNR=10*log10(Max/MSE);
disp(PSNR)
% %MSE for h2 filterd image

Im=(A-h2).^2;
i1=Im(:,:,1);
i2=Im(:,:,2);
i3=Im(:,:,3);
S1=sum(i1,'all');
S2=sum(i2,'all');
S3=sum(i3,'all');
Avg1=(S1+S2+S3)/3;
E=1/(N1*N2);
MSE1=E.*Avg1;
disp("MSE of image with h2 filter :")
disp(MSE1)

Max=255^2;
PSNR1=10*log10(Max/MSE1);
disp("PSNR for image with h2 filter :")
disp(PSNR1)



