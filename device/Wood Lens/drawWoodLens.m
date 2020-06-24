clear;
clc;

cm = 10^-2;
mm = 10^-3;
um = 10^-6;
nm = 10^-9;


grid = importdata('Grid.dat');
urd = 1.0;
erd = 1.0;

%Initialize Layers to device
UR = urd * ones(grid(4),grid(5));
ER = erd * ones(grid(4),grid(5));

save ER_Layer2.dat ER;
save UR_Layer2.dat UR;

dx = grid(2) / grid(4);
dy = grid(3) / grid(5);

R = 7 * mm;

for i = 1:grid(4)
  for j = 1:grid(5)
     r2 = ((i - grid(4)/2)*dx)^2 + ((j-grid(5)/2)*dy)^2;
     if(sqrt(r2) <= R)
        ER(i,j) = (1.5 - 5000 * r2)^2;
     end
  end
end

imagesc(ER)
colorbar;
save ER_Layer1.dat ER;
save UR_Layer1.dat UR;

% r = linspace(0,0.007,100);
% plot(r,(1.5 - 5000*(r.^2)))