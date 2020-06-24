function [R,T,m,F] = fdfd2d(lam0,UR2,ER2,RES2,NPML,kinc,pol)
% FDFD2D      Two-Dimensional Finite-Difference Frequency-Domain
%
% This MATLAB code simulates optical structures using the 
% finite-difference frequency-domain method.
%
% INPUT ARGUMENTS
% lam0 is the free space wavelength
% UR2 contains the relative permeability on a 2X grid
% ER2 contains the relative permittivity on a 2X grid
% NPML is the size of the PML on the 1X grid
% RES2 = [dx2 dy2]
% kinc is the indicent wave vector
% pol is the polarization ('E' or 'H')
%
% OUTPUT ARGUMENTS
% R contains diffraction efficiencies of reflected waves
% T contains diffraction efficiencies of transmitted waves
% m contains the indices of the harmonics in R and T
% F is the computed field
%
% Raymond C. Rumpf, Ph.D.
% Associate Professor of Electrical and Computer Engineering
% University of Texas at El Paso
% EL Paso, TX 79968
%
% Example code written for short course on
% "Introduction to Optical Simulation Using the Finite-Difference
% Frequency-Domain Method."

%% HANDLE INPUT AND OUTPUT ARGUMENTS

% DETERMINE SIZE OF GRID
[Nx2,Ny2] = size(ER2);
dx2 = RES2(1);
dy2 = RES2(2);

% 1X GRID PARAMETERS
Nx = Nx2/2;     dx = 2*dx2;
Ny = Ny2/2;     dy = 2*dy2;

% COMPUTE MATRIX SIZE
M = Nx*Ny;

% COMPUTE REFRACTIVE INDEX IN REFLECTION REGION
erref = ER2(:,1);          erref = mean(erref(:));
urref = UR2(:,1);          urref = mean(urref(:));
nref = sqrt(erref*urref);
if erref<0 && urref<0
    nref = - nref;
end

% COMPUTE REFRACTIVE INDEX IN TRANSMISSION REGION
ertrn = ER2(:,Ny2);        ertrn = mean(ertrn(:));
urtrn = UR2(:,Ny2);        urtrn = mean(urtrn(:));
ntrn = sqrt(ertrn*urtrn);
if ertrn<0 && urtrn<0
    ntrn = - ntrn;
end

%% INCORPORATE PERFECTLY MATCHED LAYER BOUNDARY CONDITION

% PML PARAMETERS
N0   = 376.73032165;                    %free space impedance
amax = 3;
cmax = 1;
p    = 3;

% INITIALIZE PML TO PROBLEM SPACE
sx = ones(Nx2,Ny2);
sy = ones(Nx2,Ny2);

% COMPUTE FREE SPACE WAVE NUMBERS
k0 = 2*pi/lam0;

% Y PML
N = 2*NPML;
for n = 1 : N
    % compute PML value
    ay = 1 + amax*(n/N)^p;
    cy = cmax*sin(0.5*pi*n/N)^2;
    s  = ay*(1-1i*cy*N0/k0);
    % incorporate value into PML
    sy(:,N-n+1)   = s;
    sy(:,Ny2-N+n) = s;
end

% COMPUTE TENSOR COMPONENTS WITH PML
ER2xx = ER2 ./ sx .* sy;
ER2yy = ER2 .* sx ./ sy;
ER2zz = ER2 .* sx .* sy;

UR2xx = UR2 ./ sx .* sy;
UR2yy = UR2 .* sx ./ sy;
UR2zz = UR2 .* sx .* sy;

% OVERLAY MATERIALS ONTO 1X GRID
ERxx = ER2xx(2:2:Nx2,1:2:Ny2);
ERyy = ER2yy(1:2:Nx2,2:2:Ny2);
ERzz = ER2zz(1:2:Nx2,1:2:Ny2);
URxx = UR2xx(1:2:Nx2,2:2:Ny2);
URyy = UR2yy(2:2:Nx2,1:2:Ny2);
URzz = UR2zz(2:2:Nx2,2:2:Ny2);

% CLEAR TEMPORARY VARIABLES
clear N0 amax cmax p sx sy n N ay cy s;
clear UR2 ER2 ER2xx ER2yy ER2zz UR2xx UR2yy UR2zz;

%% PERFORM FINITE-DIFFERENCE FREQUENCY-DOMAIN ANALYSIS

% FORM DIAGONAL MATERIAL MATRICES
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));

URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% COMPUTE DERIVATIVE OPERATORS
NS   = [Nx Ny];
RES  = [dx dy];
BC   = [-2 -2 0 0];
[DEX,DEY,DHX,DHY] = yeeder(NS,k0*RES,BC,kinc/k0);

% COMPUTE FIELD MATRIX
switch pol
    case 'E',
        A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
    case 'H',
        A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
    otherwise,
        error('Unrecognized polarization.');
end

% COMPUTE SOURCE FIELD
xa    = [0:Nx-1]*dx;
ya    = [0:Ny-1]*dy;
[Y,X] = meshgrid(ya,xa);
fsrc  = exp(-i*(kinc(1)*X+kinc(2)*Y));
fsrc  = fsrc(:);

% COMPUTE SCATTERED-FIELD MASKING MATRIX
Q = zeros(Nx,Ny);
Q(:,1:NPML+2) = 1;
Q = diag(sparse(Q(:)));

% COMPUTE SOURCE VECTOR
f = (Q*A-A*Q)*fsrc;

% PREPARE MEMORY
clear NS RES BC DEX DEZ DHX DHZ;
clear ya X Y fsrc;
clear ERxx ERyy ERzz URxx URyy URzz;

% COMPUTE FIELD
F = A\f;                  %backward division is used here!!
F = full(F);
F = reshape(F,Nx,Ny);

%% COMPUTE DIFFRACTION EFFICIENCIES

% EXTRACT REFLECTED AND TRANSMITTED WAVES
Fref = F(:,NPML+1);
Ftrn = F(:,Ny-NPML);

% REMOVE PHASE TILT
p    = exp(+1i*kinc(1)*xa');
Fref = Fref .* p;
Ftrn = Ftrn .* p;

% COMPUTE SPATIAL HARMONICS
Fref = fftshift(fft(Fref))/Nx;
Ftrn = fftshift(fft(Ftrn))/Nx;

% COMPUTE WAVE VECTOR COMPONENTS OF THE SPATIAL HARMONICS
m   = [-floor(Nx/2):floor(Nx/2)]';
kx  = kinc(1) - 2*pi*m/(Nx*dx);
kzR = conj( sqrt((k0*nref)^2 - kx.^2) );
kzT = conj( sqrt((k0*ntrn)^2 - kx.^2) );

% COMPUTE DIFFRACTION EFFICIENCY
switch pol
    case 'E',
        R =  abs(Fref).^2 .* real(kzR/kinc(2));
        T =  abs(Ftrn).^2 .* real(kzT*urref/kinc(2)/urtrn);
    case 'H',
        R =  abs(Fref).^2 .* real(kzR/kinc(2));
        T =  abs(Ftrn).^2 .* real(kzT*erref/kinc(2)/ertrn);
end

