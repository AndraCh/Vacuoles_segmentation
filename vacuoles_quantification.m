%% Andra Chincisan, Institute of Neuropathology, USZ
% 2018
%%
close all;
clc;
clear all;
%% Read image
path_dir = 'path_name'
path_save = 'save_path_name'
image_type = '.tif';
srcFiles = dir(strcat(path_dir, '*', image_type));
for i = 1 : length(srcFiles)
    image = imread (strcat(path_dir,srcFiles(i).name)); 
    imageGray = getgrayimage (image);
    imageBinary = segmentimage (imageGray);
    %figure, imshow(image(:,:,3)); title('imageVacoulesNew');
    imageVacoules = getvacoulesimage (imageBinary);
    [numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules, Area_all,Eccentricity_all] = quantifyvacoules (imageVacoules);
    imageVacoulesNew= remove_vessels (imageGray, imageVacoules);
    %figure, imshow(imageVacoulesNew); title('imageVacoulesNew');

    C = imfuse(imageGray,imageVacoulesNew);
    figure, imshow(C); title(strcat(path_dir,srcFiles(i).name));
    imwrite(mat2gray(C),strcat(path_save, srcFiles(i).name));
end
%% Get second layer of the RBG image (blue layer) as gray image
function imageGray = getgrayimage (image)
imageGray = image(:,:,2);
%figure, imshow (imageGray); title ('Gray image'); %imcontrast
end 

%% Threshold segmentation
function imageBinary = segmentimage (imageGray)
% Adjust image contrast
imageGray = imadjust(imageGray,[0.4 0.6],[]);
% Segment image
imageBinary = im2bw(imageGray,0.7);
% Remove small objects
imageBinary = bwareaopen(imageBinary,20);
%Dilation
se = strel('disk',2);
imageBinary = imdilate(imageBinary,se);
% Erosion
se = strel('disk',2);
imageBinary = imdilate(imageBinary,se);
imageBinary = imfill(imageBinary,'holes');
%figure, imshow(imageBinary);
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
%shape = [S.Perimeter].^2 ./(4*pi*[S.Area])
% Select all large objects Area> 350 pixels
%imgBinaryObjLarge = ismember(imageLabel, find(([S.Area] <=  500))); 
% Select smaller and rounder objects
imgBinaryObjSmall = ismember(imageLabel, find(([S.Area] >= 200) & ([S.Area]< 2000) &(shape>= 0.9)));% & ([S.Eccentricity]>0.5)));
% Combine images
imageVacoules = imgBinaryObjSmall;
%figure, imshow (imageVacoules); title ('Final vacoules images');
end

%% Final objects list
function [numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules,Area_all,Eccentricity_all] = quantifyvacoules (imageVacoules)
[imageLabel, numberOfVacoules] = bwlabel(imageVacoules);
 imageCC = bwconncomp(imageVacoules, 4 );
% Geometrical parameters of the image
S = regionprops(imageCC, 'Area', 'Perimeter', 'Eccentricity'); 
sumAreaVacoules = sum([S.Area]);
sumPerimeterVacoules = sum([S.Perimeter]);
Area_all = [S.Area];
Eccentricity_all= [S.Eccentricity];
end

function mat = window_image(Cx,Cy,D,imageGray)
[max_n, max_m] = size(imageGray);
D = D /2;
new_coordinates = [round(Cx-D/2),round(Cx+D/2),round(Cy-D/2),round(Cy+D/2)];

if (Cx-D)>=0 & (Cx-D)<=max_m
    new_coordinates(1) = round(Cx-D);
end
if (Cx+D)>=0 & (Cx+D)<=max_m
    new_coordinates(2) = round(Cx+D);
end
if (Cy-D)>=0 & (Cy-D)<=max_n
    new_coordinates(3) = round(Cy-D);
end
if (Cy+D)>=0 & (Cy+D)<=max_n
    new_coordinates(4) = round(Cy+D);
end
mat = zeros ( new_coordinates(4)-new_coordinates(3),new_coordinates(2)-new_coordinates(1));
for i=1:(new_coordinates(4)-new_coordinates(3))
    for j=1:(new_coordinates(2)-new_coordinates(1))
        if (new_coordinates(3)+i<= max_n) & (new_coordinates(1)+j<= max_m)
            mat(i,j) = imageGray(new_coordinates(3)+i,new_coordinates(1)+j);
        end
    end
end
end

function [imageVacoules]= remove_vessels (imageGray, imageVacoules)
S = regionprops(imageVacoules, 'centroid', 'MajorAxisLength'); 
centroids = cat(1, S.Centroid);
diameters = [S.MajorAxisLength];
[imgLabel, numberOfObject] = bwlabel(imageVacoules);

[s1, s2] = size(imageVacoules);
for noObj=1:numberOfObject
     mat = window_image(centroids(noObj,1),centroids(noObj,2),diameters(noObj),imageGray);
     [m, n] = size(mat);
     if (m>10) & (n>10)
        if mean(mean(mat))<100%mean(mean(imageGray))70
            imageVacoules(imgLabel==noObj) = 0;
        end
     end
end
[imgLabel, numberOfObject] = bwlabel(imageVacoules);
end