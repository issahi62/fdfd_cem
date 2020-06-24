clear;
clc;

grid = importdata('Grid.dat');
urd = 1.0;
erd1 = 9.0;
erd2 = 1.0;

%Initialize Layers to device
UR = urd * ones(grid(4),grid(5));
ER = erd2 * ones(grid(4),grid(5));

Nx = grid(4)/8;
Ny = grid(5)/8;
CenterX = grid(4)/2;
CenterY = grid(5)/2;

ER(CenterX-Nx:CenterX+Nx,CenterY-Ny:CenterY+Ny) = erd1;

imagesc(ER)
colorbar;
save ER_Layer1.dat ER;
save UR_Layer1.dat UR;

