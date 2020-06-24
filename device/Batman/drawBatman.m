clear;
clc;

A = imread('bat.jpg');
A = rgb2gray(A);
A = double(A);
A(A<=200) = -1000;
A(A>200) = 1;

B = 0*A + 1;
A = A';
save ER_Layer1.dat A
save UR_Layer1.dat B
