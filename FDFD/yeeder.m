function [DEX,DEY,DHX,DHY] = yeeder(grid)
% EXTRACT GRID PARAMETERS
Nx = grid.Nx;
dx = (2*pi / grid.lam0) * grid.Lx / grid.Nx;
Ny = grid.Ny;   
dy = (2*pi / grid.lam0) * grid.Ly / grid.Ny;

% DETERMINE MATRIX SIZE
M = Nx*Ny;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATRIX
DEX = sparse(M,M);

% PLACE MAIN DIAGONALS
DEX = spdiags(-ones(M,1),0,DEX);
DEX = spdiags(+ones(M,1),+1,DEX);

% CORRECT BOUNDARY TERMS (DEFAULT TO DIRICHLET)
for ny = 1 : grid.Ny-1
    neq = grid.Nx*(ny-1) + grid.Nx;
    DEX(neq,neq+1) = 0;
end

% FINISH COMPUTATION
DEX = DEX / dx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATRIX
DEY = sparse(M,M);

% PLACE MAIN DIAGONALS
DEY = spdiags(-ones(M,1),0,DEY);
DEY = spdiags(+ones(M,1),+Nx,DEY);

% FINISH COMPUTATION
DEY = DEY / dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DHX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATRIX
DHX = sparse(M,M);

% PLACE MAIN DIAGONALS
DHX = spdiags(+ones(M,1),0,DHX);
DHX = spdiags(-ones(M,1),-1,DHX);

% FINISH COMPUTATION
DHX = DHX / dx;
for ny = 2 : grid.Ny
    neq = grid.Nx*(ny-1) + 1;
    DHX(neq,neq-1) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DHY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATRIX
DHY = sparse(M,M);

% PLACE MAIN DIAGONALS
DHY = spdiags(+ones(M,1),0,DHY);
DHY = spdiags(-ones(M,1),-Nx,DHY);

% FINISH COMPUTATION
DHY = DHY / dy;

