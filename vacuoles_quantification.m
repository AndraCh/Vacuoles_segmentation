%% Andra Chincisan, Institute of Neuropathology, USZ
% 2018
%%
close all;
clc;
clear all;
%% Read image
name = 'test_image.tif';
image = imread (name);
imageGray = getgrayimage (image);
imageBinary = segmentimage (imageGray);
imageVacoules = getvacoulesimage (imageBinary);
[numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules] = quantifyvacoules (imageVacoules)

%% Get second layer of the RBG image (blue layer) as gray image
function imageGray = getgrayimage (image)
imageGray = image(:,:,2);
figure, imshow (imageGray); title ('Gray image'); 
end 

%% Threshold segmentation
function imageBinary = segmentimage (imageGray)
% Adjust image contrast
imageGray = imadjust(imageGray,[0.5 0.9],[]);
% Segment image
imageBinary = im2bw(imageGray,0.8); 
% Remove small objects
imageBinary = bwareaopen(imageBinary,10);
% Dilation
se = strel('disk',2);
imageBinary = imdilate(imageBinary,se);
% Erosion
se = strel('disk',2);
imageBinary = imdilate(imageBinary,se);
figure, imshow(imageBinary);
end

%% Get vacoules
function imageVacoules = getvacoulesimage (imageBinary)
%  Connected components
imageCC = bwconncomp(imageBinary, 4 );
% Geometrical parameters of the image
S = regionprops(imageCC, 'Area', 'Perimeter', 'Eccentricity', 'EquivDiameter'); 
% Label objects in the image
imageLabel = labelmatrix(imageCC);
% Calculate object's shape
shape = 4*pi*[S.Area]./power ([S.Perimeter],2);
% Select all large objects Area> 350 pixels
imgBinaryObjLarge = ismember(imageLabel, find(([S.Area] >=  350) &(shape>= 0.5))); 
% Select smaller and rounder objects
imgBinaryObjSmall = ismember(imageLabel, find(([S.Area] >= 200) & ([S.Area]< 350) &(shape>= 0.9) & ([S.Eccentricity]<=0.9)));
% Combine images
imageVacoules = imadd(im2double(imgBinaryObjLarge),im2double(imgBinaryObjSmall));
figure, imshow (imageVacoules); title ('Final vacoules images');
end

%% Final objects list
function [numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules] = quantifyvacoules (imageVacoules)
[imageLabel, numberOfVacoules] = bwlabel(imageVacoules);
 imageCC = bwconncomp(imageVacoules, 4 );
% Geometrical parameters of the image
S = regionprops(imageCC, 'Area', 'Perimeter'); 
sumAreaVacoules = sum([S.Area]);
sumPerimeterVacoules = sum([S.Perimeter]);
end
