clear;
clc;

%Main Function for FDFD
%%  Simulation Parameters
fn = '../device/Free Space/';
pmlFlag = true;
pol = 'TE';

%%  Grid Calculation
grid = setupGrid(strcat(fn,'Grid.dat'));
NLambda = grid.lam0 / grid.dx;

%%  Device Calculation (UPML)
 percentPML = 10 * 0.01;
 pmlX = round(percentPML * grid.Nx);
 pmlY = round(percentPML * grid.Ny);

[pml.sx, pml.sy] = calcpml(grid,pmlX,pmlY);

if(~pmlFlag)
    pml.sx = 1;
    pml.sy = 1;
end

device = [];
device = setupDevice(fn,device,pml);

%%  Derivative Operator
%Derivative Matrices
[DEX,DEY,DHX,DHY] = yeeder(grid);

switch pol
    case 'TM'
        A = DHX/device.URyy*DEX + DHY/device.URxx*DEY + device.ERzz;
    case 'TE'
        A = DEX/device.ERyy*DHX + DEY/device.ERxx*DHY + device.URzz;
    otherwise
        error('Unrecognized polarization.');
end
%%  Source Calculation
src = setupSrc(grid,A,pmlX,pmlY);

%%  Compute Field
Psi = A\src;                  
Psi = full(Psi);
Psi = reshape(Psi,grid.Nx,grid.Ny);

Pref = Psi(:,pmlY+1);
Ptrn = Psi(:,grid.Ny - pmlY);
phase =  exp(+1i * (2*pi / grid.lam0) * [0:grid.Nx]'*grid.dx);
%Pref = Pref .* phase;
%Ptrn = Ptrn .* phase;


imagesc([0:grid.Nx-1]*grid.dx,[0:grid.Ny-1]*grid.dy,real(Psi)');
xlabel('x (\mum)');
ylabel('y (\mum)');
title('E_z FIELD');
colorbar;

