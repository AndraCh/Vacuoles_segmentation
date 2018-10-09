# Vacuoles_segmentation
Matlab functions to segment and quantify vacuoles in histological images.

The algorithm was developed in MATLAB R2016B, using the Image Processing toolbox. 
Image segmentation was performed using Otsu’s thresholding method (T=0.7). 
Only round vacuoles with an area in the range of [200 pixels, 2000 pixels] and with a shape values > 0.9 
(where shape was calculated as shape = (4*PI*Area) / (Perimeter^2)) were quantified. 


