%% Andra Chincisan, Institute of Neuropathology, USZ
% 2018
%%
close all;
clc;
clear all;
%% Read image
path_dir = 'path_dir\'
path_save = 'path_save\'
path_xls = 'path_save'
save_xls = 'file.xls'
image_type = '.tif';
srcFiles = dir(strcat(path_dir, '*', image_type));
for i = 1 : length(srcFiles)
    close all;
    % Read image
    image = imread (strcat(path_dir,srcFiles(i).name)); 
    % Gray image
    imageGray = getgrayimage (image);
    % Segment image
    imageBinary = segmentimage (imageGray);
    figure, imshow(image(:,:,3)); title('imageVacoulesNew');
    % Get vacuoles
    imageVacoules = getvacoulesimage (imageBinary);
    % Process vacuoles
    imageVacoulesNew= remove_vessels (imageGray, imageVacoules);
    % Quantify vacuoles
    [numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules, Area_all,Eccentricity_all] = quantifyvacoules (imageVacoulesNew);
    figure, imshow(imageVacoulesNew); title('imageVacoulesNew');
    % Fuse original and segmented image
    C = imfuse(imageGray,imageVacoulesNew);
    figure, imshow(C); title(strcat(path_dir,srcFiles(i).name));
    % Save fused image
    imwrite(mat2gray(C),strcat(path_save, srcFiles(i).name));
    name_image =string (srcFiles(i).name);
    % Create arry with data for each each
    all_param(i, :) = [{char(name_image)} numberOfVacoules sumAreaVacoules sumPerimeterVacoules Area_all]
end
% Save data in Excel file
xlswrite(strcat (path_xls, '\het_cortex'), save_xls)

%% Get second layer of the RBG image (blue layer) as gray image
function imageGray = getgrayimage (image)
imageGray = image(:,:,2);
%figure, imshow (imageGray); title ('Gray image'); %imcontrast
end 

%% Threshold segmentation
function imageBinary = segmentimage (imageGray)
% Adjust image contrast
imageGray = imadjust(imageGray,[0.2 0.8],[]);  
% Binarize image
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
% Select smaller and rounder objects
imgBinaryObjSmall = ismember(imageLabel, find(([S.Area] >= 200) & ([S.Area]< 2000) &(shape>= 0.9)));% & ([S.Eccentricity]>0.5)));
% Combine images
imageVacoules = imgBinaryObjSmall;
%figure, imshow (imageVacoules); title ('Final vacoules images');
end

%% Final objects list
function [numberOfVacoules,sumAreaVacoules,sumPerimeterVacoules,mean_area_all,mean_eccentricity_all] = quantifyvacoules (imageVacoules)
[imageLabel, numberOfVacoules] = bwlabel(imageVacoules);
 imageCC = bwconncomp(imageVacoules, 4 );
% Geometrical parameters of the image
S = regionprops(imageCC, 'Area', 'Perimeter', 'Eccentricity'); 
% Calculate sum of all vacuole areas in an image
sumAreaVacoules = sum([S.Area]);
% Calculate sum of all vacuole perimeters in an image
sumPerimeterVacoules = sum([S.Perimeter]);
% Calculate mean area of vacuole in an image
mean_area_all = mean([S.Area]);
% Calculate mean eccentricity of vacuole in an image
mean_eccentricity_all= mean([S.Eccentricity]);
end

function mat = window_image(Cx,Cy,D,imageGray)
[max_n, max_m] = size(imageGray);
D = D/3;
new_coordinates = [round(Cx-D/2),round(Cx+D/2),round(Cy-D/2),round(Cy+D/2)];
for i=1:length(new_coordinates)
    if new_coordinates(i) <=0
        new_coordinates(i) = 0;  
    end
end
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
% Create a small image for each of the initial objects that are considered
% vacuoles. If dark regions appear in the close proximity of an a object it 
% means that there is a vessel close to the object and that this object is 
% not a real vacuole, so it should be remove from the vacuoles list.
% Compare the average intensity of small images with a threshold T =
% 150. Objects with an average intensity <150 will be removed.
 
% Image parameters
S = regionprops(imageVacoules, 'centroid', 'MajorAxisLength'); 
centroids = cat(1, S.Centroid);
diameters = [S.MajorAxisLength];
[imgLabel, numberOfObject] = bwlabel(imageVacoules);

[s1, s2] = size(imageVacoules);
for noObj=1:numberOfObject
     mat = window_image(centroids(noObj,1),centroids(noObj,2),diameters(noObj),imageGray);
     [m, n] = size(mat);
     if (m>10) & (n>10)
        if mean(mean(mat))<150
            imageVacoules(imgLabel==noObj) = 0;
        end
     end
end
% Count objects
[imgLabel, numberOfObject] = bwlabel(imageVacoules);
end