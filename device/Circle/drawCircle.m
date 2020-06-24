clear;
clc;

grid = importdata('Grid.dat');
urd = 1.0;
erd1 = -1000;
erd2 = 1.0;

%Initialize Layers to device
UR = urd * ones(grid(4),grid(5));
ER = erd2 * ones(grid(4),grid(5));

[cols rows] = meshgrid(1:grid(4), 1:grid(5));
% Next create the circle in the image.
centerX = round(grid(4)/2);
centerY = round(grid(5)/2);
radius = round(grid(4)/4);
circlePixels = (rows - centerY).^2 ...
    + (cols - centerX).^2 <= radius.^2;
ER(circlePixels) = erd1;

imagesc(ER)

save ER_Layer1.dat ER;
save UR_Layer1.dat UR;

