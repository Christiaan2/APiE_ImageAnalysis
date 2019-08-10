%{
    Name: Christiaan Reurslag
    Studentnummer: S1495089
    Assignment: Image Analysis Assignment 2
    MATLAB R2016b
%}

%Close all figures and clear workspace and command window
clear all
close all
clc

%% 1a
%Load images into workspace
pattern1 = imread('pattern1.tif');
pattern2 = imread('pattern2.tif');

%Calculate and plot cross-correlation of pattern1 and itself
Cross1 = myCrossCorrelation(pattern1,pattern1);
figure
surf(Cross1);

%Calculate and plot cross-correlation of pattern2 and itself
Cross2 = myCrossCorrelation(pattern2,pattern2);
figure
surf(Cross2);

%% 1c
%Calculate and plot cross-correlation of pattern1 and pattern2
Cross12 = myCrossCorrelation(pattern1,pattern2);
figure
surf(Cross12);

%% 1d
%Calculate the displacment
[i,j] = find(Cross12 == max(max(Cross12)));
Displacement = ceil(size(Cross12)/2)-[i,j];

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