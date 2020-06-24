clear;
clc;

A = imread('superman.png');
A = double(A);
A(A<=200) = -1000;
A(A>200) = 1;

A = A';

aMat = ones(1024,1024);
aMat(200:200+600-1,200:200+460-1) = A;

A = aMat(1:4:end,1:4:end);

B = 0*A + 1;

save ER_Layer1.dat A
save UR_Layer1.dat B
