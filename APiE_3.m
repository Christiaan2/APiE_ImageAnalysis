%{
    Name: Christiaan Reurslag
    Studentnummer: S1495089
    Assignment: Image Analysis Assignment 3
    MATLAB R2016b
%}

%Close all figures and clear workspace and command window
clear all
close all
clc

%% 1a
%Load images into workspace
E1 = imread('E001_1.tif');
E2 = imread('E001_2.tif');

%Perform adaptive histogram equalisation
M = 16; %Size window for adaptive histogram equalisation (horizontal)
N = M;  %(vertical)
bin = linspace(0,255,256);

%Plot cdf before and after histogram equalisation
%Image E001_1.tif
E1_EQ = myAdaptiveHistEQ(E1,bin,M,N);%Adaptive histrogram equalisation
p1 = myHistogram(E1,bin);%Histogram of original image
cdf1 = myCalculateCDF(p1);%CDF of original image
p1_EQ = myHistogram(E1_EQ,bin);%Histogram of image after equalisation
cdf1_EQ = myCalculateCDF(p1_EQ);%CDF of image after equalisation
figure
hold on
plot([bin(1),repelem(bin(2:end),2),bin(end)],repelem(cdf1,2))
plot([bin(1),repelem(bin(2:end),2),bin(end)],repelem(cdf1_EQ,2))
title('CDF of image (E1)')
ylabel('CDF')
xlabel('Intensity')
legend('Before equalisation','After equalisation')
grid on
grid minor

%Image E001_2.tif
E2_EQ = myAdaptiveHistEQ(E2,bin,M,N);
p2 = myHistogram(E2,bin);
cdf2 = myCalculateCDF(p2);
p2_EQ = myHistogram(E2_EQ,bin);
cdf2_EQ = myCalculateCDF(p2_EQ);
figure
hold on
plot([bin(1),repelem(bin(2:end),2),bin(end)],repelem(cdf2,2))
plot([bin(1),repelem(bin(2:end),2),bin(end)],repelem(cdf2_EQ,2))
title('CDF of image (E2)')
ylabel('CDF')
xlabel('Intensity')
legend('Before equalisation','After equalisation')
grid on
grid minor

%Show original image and image after equalisation
figure
imshow(E1)
title('E1 before equalisation')
figure
imshow(E1_EQ/(length(bin)-1))
title('E1 after equalisation')
figure
imshow(E2)
title('E2 before equalisation')
figure
imshow(E2_EQ/(length(bin)-1))
title('E2 after equalisation')

%Calculate velocity field
N2 = 16;    %Size of sub-images (horizontal as well as vertical)
v = myVelocityField(E1_EQ,E2_EQ,N2);

%Plot velocity field
figure
ii = 1:N:length(E1_EQ(:,1));
jj = 1:N:length(E1_EQ(1,:));
Speed = NaN(size(E1_EQ));
Speed(v(1,1,3):v(end,1,3),v(1,1,4):v(1,end,4)) = sqrt(v(:,:,2).^2+v(:,:,1).^2);
SpeedPlot = imagesc(Speed);
set(SpeedPlot,'AlphaData',~isnan(Speed))
colorbar
axis equal
hold on
quiver(v(ii,jj,4),v(ii,jj,3),v(ii,jj,2),v(ii,jj,1),'r')
title('Velocity field')
xlabel('Pixel number')
ylabel('Pixel number')

function v = myVelocityField(frame1,frame2,N)
%%Calculate velocity field using particle imaging velocimetry
%   Inputs: frame1 = image at time t = t0
%           frame2 = image at time t = t0+dt
%           N = size of sub-image(horizontal and vertical)

%Preallocate the matrix v to store the displacement and location
v = zeros([size(frame1)-16+1,4]);
ii = 1:113;
jj = 1:241;
for i = 1:length(ii)
    for j = 1:length(jj)
        Cross = myCrossCorrelation(frame1(ii(i):ii(i)+N-1,jj(j):jj(j)+N-1),frame2(ii(i):ii(i)+N-1,jj(j):jj(j)+N-1));
        [imax,jmax] = find(Cross == max(max(Cross)));%find maximum of cross-correlation
        v(i,j,1:2) = ceil(size(Cross)/2)-[imax,jmax];%displacement
        v(i,j,3:4) = [ii(i)-1+N/2,jj(j)-1+N/2];%location of displacement
    end
end

end

function [p] = myHistogram(X,bin)
%%Returns the probability density function(pdf)
%   Inputs: X = matrix to transform into pdf
%           bin = vector which specifies the edges of the bins

%Preallocate the vector p
p = zeros(1,length(bin));
for i = 1:length(bin)-1
    A = bin(i)*ones(size(X));
    B = bin(i+1)*ones(size(X));
    p(i) = length(find(A <= X & X < B))/numel(X);
    if(i == length(bin)-1)
        %Only the last bin includes the edge
        p(i) = p(i) + length(find(X == B))/numel(X);
    end
end
end

function cdf = myCalculateCDF(p)
%%Calculates the cumulative distribution function
%Inputs:    p = probability density function
cdf = zeros(size(p));

for i = 1:length(cdf)
    cdf(i) = sum(p(1:i));
end
end

function I = myAdaptiveHistEQ(image,bin,M,N)
%%Perform adaptive histrogram equalisation using sliding window approach
%Inputs:    image = matrix corresponding to image
%           bin = vector which specifies the edges of the bins
%           M = Size of window (vertical)
%           N = Size of window (horizontal)

%Preallocate the matrix I to store the image after equalisation
I = zeros(size(image));
for i = 1:length(image(:,1))-M
    for j = 1:length(image(1,:))-N
        p = myHistogram(image(i:i+N,j:j+M),bin);
        cdf = myCalculateCDF(p);
        [~,I(i:i+N,j:j+M)] = myHistEQ(image(i:i+N,j:j+M),bin,cdf);
    end
end
end

function [h,I] = myHistEQ(image,bin,cdf)
%%Perform histogram equlaistion 
%Inputs:    image = image (or part of image) to perform equalisation on
%           bin = vector which specifies the edges of the bins
%           cdf = cumulative distribution function
h = zeros(size(cdf));
cdfmin = min(cdf(cdf > 0));%find non-zero minimum of cdf
%Determine new intensity values
for i = 1:length(cdf)
    h(i) = (cdf(i)-cdfmin)/(cdf(end)-cdfmin)*(length(cdf)-1);
end

%Preallocate output array containing new image
I = zeros(size(image));
for i = 1:length(bin)-1
    A = bin(i)*ones(size(image));
    B = bin(i+1)*ones(size(image));
    I(A <= image & image < B) = h(i);
    if(i == length(bin)-1)
        %Only the last bin includes the edge
        I(image == B) = h(i);
    end
end
end

function CrossCorrelation = myCrossCorrelation(image,mask)
%%Calculate the cross-correlation
%   Input: image = matrix corresponding to image
%          mask = matrix corresponding to mask

%Preallocate the matrix CrossCorrelation
CrossCorrelation = zeros(size(image)+size(mask)-1);
%Transform input matrices to double (otherwise maximum value of product is
%255)
image = double(image);
mask = double(mask);

%Step through matrix and calculate un-normalized correlation of the overlap
%between the matrices image and mask
for i = 1:length(CrossCorrelation(:,1))
    for j = 1:length(CrossCorrelation(1,:))
        maxXimage = min(j,length(image(1,:)));
        minXimage = max(j-length(mask(1,:))+1,1);
        maxYimage = min(i,length(image(:,1)));
        minYimage = max(i-length(mask(:,1))+1,1);
        maxXmask = min(length(mask(1,:)),length(CrossCorrelation(1,:))-j+1);
        minXmask = max(length(mask(1,:))-j+1,1);
        maxYmask = min(length(mask(:,1)),length(CrossCorrelation(:,1))-i+1);
        minYmask = max(length(mask(:,1))-i+1,1);
        CrossCorrelation(i,j) = sum(sum((image(minYimage:maxYimage,minXimage:maxXimage) - ...
            mean2(image(minYimage:maxYimage,minXimage:maxXimage))).* ...
            (mask(minYmask:maxYmask,minXmask:maxXmask) - ...
            mean2(mask(minYmask:maxYmask,minXmask:maxXmask)))));
    end
end
end