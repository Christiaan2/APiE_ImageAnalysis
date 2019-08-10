%{
    Name: Christiaan Reurslag
    Studentnummer: S1495089
    Assignment: Image Analysis Assignment 1
    MATLAB R2016b
%}

%Close all figures and clear workspace and command window
clear all
close all
clc

%% 1a
%Load image into workspace
RGB = imread('bladcellen.png');

%Transform colour image to grayscale image
grayScale = myRGB2Gray(RGB);

%Plot histrogram to check if we have indeed a bimodal distribution
bin = linspace(0,255,256);
p = myHistogram(grayScale,bin);
figure
area([bin(1),repelem(bin(2:end),2),bin(end)],repelem(p,2))
title('Histogram of image in grayscale')
ylabel('PDF')
xlabel('Intensity')
grid on
grid minor

%Separate cells from background using Otsu's method
[binary, threshold] = myOTSU(grayScale,bin,p);
disp(['Treshold value from OTSU''s method: ', num2str(threshold)])
%Show binary image
figure
imshow(binary)

%Calculate number of cells using two pass filter
[TwoPass,N] = myTwoPassAlgorithm(binary,25);
%show image with every cell having it's own colour
figure
imagesc(TwoPass)
axis equal
disp(['Number of cells: ', num2str(N)])

%% 1b
%Determine cell numbers which don't touch the border of the image
CellsnotTouching = myCellsnotTouching(TwoPass,N);

%Determine area of cells in CellsnotTouching
Area = myAreaCells(TwoPass,CellsnotTouching);

%Determine center of mass of cells in CellsnotTouching
CenterofMass = myCenterofMass(TwoPass,Area);

%Determine major and minor axes of cells in CellsnotTouching
[labda,theta] = myMajorMinorAxes(TwoPass,Area,CenterofMass);

%Determine eccentricity of cells in CellsnotTouching
epsilon = myEccentricity(labda);

%Show original images with ellipses and label
figure
imshow(RGB)
hold on
plot(CenterofMass(:,2),CenterofMass(:,3),'r+','MarkerSize',8,'LineWidth',1.5)
text(CenterofMass(:,2)+1,CenterofMass(:,3)-4,num2str(CenterofMass(:,1)),'Color','k','FontSize',14)
for i = 1:length(CellsnotTouching)
    [x,y] = myEllipse(CenterofMass(i,2),CenterofMass(i,3),2*sqrt(labda(i,3)),2*sqrt(labda(i,2)),theta(i,2),500);
    plot(x,y,'r','LineWidth',1.5)
end


function grayScale = myRGB2Gray(RGB)
%%Transform RGB image to grayScale image (same as Matlab rgb2gray function)
%   Inputs: RGB -> RGB(:,:,1) = Red intensity, RGB(:,:,2) = Green intensity,
%                  RGB(:,:,3) = Blue intensity

%Transform uint8 datatype to double
RGB = double(RGB);

%calculate grayscale
RGB(:,:,1) = RGB(:,:,1) * 0.2989;
RGB(:,:,2) = RGB(:,:,2) * 0.5870;
RGB(:,:,3) = RGB(:,:,3) * 0.1140;
grayScale = RGB(:,:,1)+RGB(:,:,2)+RGB(:,:,3);
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

function [binary, threshold] = myOTSU(grayScale,bin,p)
%%Determine threshold value using OTSU's method and transform image to
%%binary image
%   Inputs: grayscale = matrix used by OTSU's method
%           bin = vector which specifies the edges of the bins
%           p = probability density function

%Preallocate array to store intra-class variance
sigma2_w = zeros(1,length(bin)-1);

for i = 1:length(bin)-1
    w_F = sum(p(1:i)); %Class probability foreground
    u_F = sum(p(1:i).*bin(1:i))/w_F; %Mean intensity forground
    sigma2_F = sum(p(1:i).*(bin(1:i)-u_F).^2); %Variance of foreground
    w_B = sum(p(i+1:end)); %Class probability background
    u_B = sum(p(i+1:end).*bin(i+1:end))/w_B; %Mean intensity background
    sigma2_B = sum(p(i+1:end).*(bin(i+1:end)-u_B).^2); %Variance of background
    sigma2_w(i) = w_F.*sigma2_F+w_B.*sigma2_B;
end

%Determine bin value (threshold) where intra-class variance has a minimum
minimum = min(sigma2_w);
threshold = bin(sigma2_w == minimum);
%Plot intra-class variance to verify minimum value
if(1)    
    figure
    plot(linspace(bin(1),bin(end-1),length(bin)-1),sigma2_w)
    xlabel('Intensity')
    ylabel('intra-class variance')
    title('Determine mininum to find threshold value')
    grid on
    grid minor
end
%Transform greyscale to binary image
binary = false(size(grayScale));
binary(grayScale <= threshold) = true;
binary(grayScale > threshold) = false;
end

function [TwoPass,N] = myTwoPassAlgorithm(binary,minSize)
%%Calculate number of foreground elements using two pass algorithm
%   Inputs: binary = binary image
%           minSize = minimum size (pixels) of foreground elements

%Preallocate TwoPass matrix with NaNs
TwoPass = NaN(size(binary));
count = 0; %Keeps track of next unused number

%First pass
for i = 1:length(binary(:,1))
    for j = 1:length(binary(1,:))
        if(binary(i,j))
            l = myFindSmallestNeighbour(TwoPass,i,j);
            if(isnan(l))
                count = count + 1;
                l = count;
            end
            TwoPass(i,j) = l;
        end
    end
end
k=0; %keeps track of number of iterations
TwoPassOld = TwoPass;   % When TwoPassOld = TwoPass, iterative process is stopped
%Second pass (I found out that one iteration for the second pass is not
%enough)
while(true)
    for i = length(binary(:,1)):-1:1
        for j = length(binary(1,:)):-1:1
            if(~isnan(TwoPass(i,j)))
                TwoPass(i,j) = myFindSmallestNeighbour(TwoPass,i,j);
            end
        end
    end
    k = k+1;
    if(sum(sum(TwoPassOld == TwoPass)) == sum(sum(binary)))
        %When TwoPassOld = TwoPass, iterative process is stopped
        disp(['Number of iterations: ', num2str(k)])
        break;
    end
    TwoPassOld = TwoPass;
end
TwoPass(isnan(TwoPass)) = 0; %Change all NaNs to 0
NCell = unique(TwoPass); %Count number of foreground elements
NCell = NCell(2:end); %First 0 element is not a foreground element

%Remove foreground elements smaller than minSize and make nummering
%continous
N = 0; %Number of foreground elements
for i = 1:length(NCell)
    if(sum(sum(TwoPass == NCell(i))) < minSize)
        TwoPass(TwoPass == NCell(i)) = 0;
    else
        N = N+1;
        TwoPass(TwoPass == NCell(i)) = N;
    end
end
end

function m = myFindSmallestNeighbour(M,i,j)
%%Returns the smallest neighbour of element (i,j in matrix M
%   Inputs: M = Matrix
%           i = row number
%           j = column number

%coordinates of all nine neighbours
c = [i-1,j-1; i-1,j; i-1,j+1;
    i,j-1; i,j; i,j+1;
    i+1,j-1; i+1,j; i+1,j+1];
m = NaN;
%Loop through all neighbours and return smallest element
for k = 1:length(c(:,1))
    if(c(k,1) < 1 ||c(k,2) < 1 || c(k,1) > length(M(:,1)) || c(k,2) > length(M(1,:)))
        continue
    end
    m = min(M(c(k,1),c(k,2)),m);
end
end

function CellsnotTouching = myCellsnotTouching(TwoPass,N)
%%Determine cell numbers which are not touching the border of the image
%   Inputs: TwoPass = Matrix returned by myTwoPassAlgorithm
%           N = Number of foreground elements

%Preallocate array
CellsnotTouching = linspace(1,N,N);

%When pixel of foreground element is at the side of the image, number of
%foreground element is removed from array CellsnotTouching
for i = 1:length(TwoPass(:,1))
    if(ismember(TwoPass(i,1),CellsnotTouching))
        CellsnotTouching(CellsnotTouching == TwoPass(i,1)) = [];
    end
    if(ismember(TwoPass(i,end),CellsnotTouching))
        CellsnotTouching(CellsnotTouching == TwoPass(i,end)) = [];
    end
end
for i = 1:length(TwoPass(1,:))
    if(ismember(TwoPass(1,i),CellsnotTouching))
        CellsnotTouching(CellsnotTouching == TwoPass(1,i)) = [];
    end
    if(ismember(TwoPass(end,i),CellsnotTouching))
        CellsnotTouching(CellsnotTouching == TwoPass(end,i)) = [];
    end
end
end

function Area = myAreaCells(TwoPass,CellsnotTouching)
%%returns the Area (number of pixels) of foreground elements in
%%CellsnotTouching
%   Inputs: TwoPass = Matrix returned by myTwoPassAlgorithm
%           CellsnotTouching = array containing the numbers of
%           foreground elements

%Preallocate matrix Area
Area = zeros(length(CellsnotTouching),2);
Area(:,1) = CellsnotTouching;
for i=1:length(CellsnotTouching)
    Area(i,2) = sum(sum(TwoPass == CellsnotTouching(i)));
end
end

function CenterofMass = myCenterofMass(TwoPass,Area)
%%Returns the center of mass of foreground elements
%   Inputs: TwoPass = Matrix returned by myTwoPassAlgorithm
%           Area = Matrix containing number foreground element and area of
%           foreground element

%Preallocate CenterofMass matrix
CenterofMass = zeros(length(Area(:,1)),3);
CenterofMass(:,1) = Area(:,1);
%Determine Center of mass
for i=1:length(Area(:,1))
    [I,J] = find(TwoPass == Area(i,1));
    CenterofMass(i,2) = sum(J)/Area(i,2); %horizontal direction
    CenterofMass(i,3) = sum(I)/Area(i,2); %vertical direction
end
end

function [labda,theta] = myMajorMinorAxes(TwoPass,Area,CenterofMass)
%%Determine major axis, minor axis and angle of ellipse surroundig
%%foreground elements
%   Inputs: TwoPass = Matrix returned by myTwoPassAlgorithm
%           Area = Matrix containing number foreground element and area of
%           foreground element
%           CenterofMass = Matrix containing number foreground element and
%           coordinates center of mass

%Preallocate arrays
labda = zeros(length(Area(:,1)),3);
labda(:,1) = Area(:,1);
theta = zeros(length(Area(:,1)),2);
theta(:,1) = Area(:,1);

for i = 1:1:length(Area(:,1))
    [I,J] = find(TwoPass == Area(i,1));
    %second order moments to create covariance matrix
    u_20 = sum(J.^2)/Area(i,2)-CenterofMass(i,2)^2;
    u_02 = sum(I.^2)/Area(i,2)-CenterofMass(i,3)^2;
    u_11 = sum(I.*J)/Area(i,2)-CenterofMass(i,2)*CenterofMass(i,3);
    %eigenvalues of covariance matrix
    labda(i,2) = 0.5*((u_20+u_02)+sqrt(4*u_11^2+(u_20-u_02)^2));
    labda(i,3) = 0.5*((u_20+u_02)-sqrt(4*u_11^2+(u_20-u_02)^2));
    %orientation
    theta(i,2) = 0.5*atan(2*u_11/(u_20-u_02));
end
end

function epsilon = myEccentricity(labda)
%%Return eccentricity
%   Inputs: labda = eigenvectors of covariance matrix
epsilon(:,1) = labda(:,1);
epsilon(:,2) = sqrt(1-labda(:,3)./labda(:,2));
end

function [x,y] = myEllipse(x0,y0,a,b,theta,N)
%%return x and y coordinates used to drawn an ellipse
%   inputs: x0 = center of ellipse (horizontal)
%           y0 = center of ellipse (vertical
%           a = minor axis
%           b = major axis
%           theta = orientation of ellipse
%           N = number of coordinates to approach ellipse

t = linspace(0,1,N);
x = a*cos(2*pi*t)*cos(theta)-b*sin(2*pi*t)*sin(theta)+x0;
y = a*cos(2*pi*t)*sin(theta)+b*sin(2*pi*t)*cos(theta)+y0;
end