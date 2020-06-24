clear;
clc;

%Main Function for FDFD
%%  Simulation Parameters
fn = '../device/Free space/';
pmlFlag = true;
srcAngle = 0;
pol = 'TM';

%%  Grid Calculation
grid = setupGrid(strcat(fn,'Grid.dat'));
NLambda = grid.lam0 / grid.dx;

%%  Device Calculation (CPML)
 percentPML = 10 * 0.01;
 pmlX = round(percentPML * grid.Nx);
 pmlY = round(percentPML * grid.Ny);

[pml.sx, pml.sy] = calcpml(grid,pmlX,pmlY);

if(~pmlFlag)
    pml.sx = 1;
    pml.sy = 1;
end

device = [];
tic
device = setupDevice(fn,device,pml);

Sx = diag(sparse(1 ./ pml.sx(:)));
Sy = diag(sparse(1 ./ pml.sy(:)));

%%  Derivative Operator
%Derivative Matrices
[DEX,DEY,DHX,DHY] = yeeder(grid);

switch pol
    case 'TM'
        A = Sx*DHX/device.URyy*Sx*DEX + Sy*DHY/device.URxx*Sy*DEY + device.ERzz;       
    case 'TE'
        A = Sx*DEX/device.ERyy*Sx*DHX + Sy*DEY/device.ERxx*Sy*DHY + device.URzz;
    otherwise
        error('Unrecognized polarization.');
end

%%  Source Calculation
src = setupSrc(grid,A,pmlX,pmlY,srcAngle);

%%  Compute Fields
Psi = A\src.vec;
toc
Psi = full(Psi);
Psi = reshape(Psi,grid.Nx,grid.Ny);

%Total Field
%Psi = Psi + src.SF;

subplot(2,2,3)
imagesc([0:grid.Nx-1]*grid.dx / (10e-6),[0:grid.Ny-1]*grid.dy / (10e-6),real(Psi)');
xlabel('x (\mum)');
ylabel('y (\mum)');
title('E Field (V/m)');
colorbar;

subplot(2,2,4)
imagesc([0:grid.Nx-1]*grid.dx / (10e-6),[0:grid.Ny-1]*grid.dy / (10e-6),(abs(Psi)).^2');
xlabel('x (\mum)');
ylabel('y (\mum)');
title('Intensity W / m^2');
colorbar;

%% COMPUTE Scattered Fields
% EXTRACT REFLECTED / TRANSMITTED WAVES
k0 = 2*pi / grid.lam0;
refX = pmlX+1;
refY = pmlY+1;
trnX = grid.Nx-pmlX;
trnY = grid.Ny-pmlY;

PrefX = (abs(Psi(refX,refY:trnY)).^2)*grid.dy;
PrefY = (abs(Psi(refX:trnX,refY)).^2)*grid.dx;
PtrnX = (abs(Psi(trnX,refY:trnY)).^2)*grid.dy;
PtrnY = (abs(Psi(refX:trnX,trnY)).^2)*grid.dx;

PscatdB = 10*log10(sum(PrefX) + sum(PrefY) + sum(PtrnX) + sum(PtrnY));

