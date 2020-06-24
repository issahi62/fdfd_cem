clear;
clc;
%clf;
%%
%   Cylindrical Test
%
    freq = 2 * 10^9;
    cc = 3*10^8;
    lam0 = cc / freq;
    
    k_rho = 2*pi / lam0;
    Nlambda = 20;
   % Nlambda = 2;
    sizeLambda = 5; 
   sizeR = sizeLambda * lam0;
    numR = Nlambda * sizeLambda;
    d_rho = sizeR / numR;
    rho = d_rho:d_rho:sizeR;
    invRho = 1 ./ rho;
    
    grid.Nx = numR;
    grid.Lx = sizeR;
    grid.Ny = round(grid.Nx);% - 0.3*grid.Nx);
    grid.Ly = grid.Lx;
    grid.lam0 = lam0;
    
    [DeZ,DeR,DhZ,DhR] = yeeder(grid);
    %invRho = meshgrid(invRho)';
    %invRho = diag(sparse(invRho(:)));
    invRhoTmp = sparse(eye(grid.Nx*grid.Ny));
    jj = 1;
    s = 0;
    for ii = 1:grid.Nx*grid.Ny
       invRhoTmp(ii,ii) = invRho(jj);
       if(jj == length(invRho))
           jj = 0;
       end
       jj = jj+1;
    end
    invRho = invRhoTmp;
    
    material = ones(grid.Nx,grid.Ny);
    urr = diag(sparse(material(:)));
    uzz = urr;
    epp = urr;
    
    %PML
    percentPML = 10 * 0.01;
    pmlX = round(percentPML * grid.Nx);
    pmlY = round(percentPML * grid.Ny);

    [pml.sx, pml.sy] = calcpml(grid,pmlX,pmlY);
    Sx = diag(sparse(1 ./ pml.sx(:)));
    Sy = diag(sparse(1 ./ pml.sy(:)));
    
    %Source
    za    = [0:grid.Nx-1]*d_rho;
    ra    = [0:grid.Ny-1]*d_rho;
    [Z,R] = meshgrid(za,ra);
    fsrc  = exp(-1i*(k_rho*R));
       
    % COMPUTE SCATTERED-FIELD MASKING MATRIX
    Q = zeros(grid.Nx,grid.Ny);
    Q(:,1:pmlY+2) = 1; %R
    %Q(1:pmlX+2,:) = 1;  %Z
    Q(:,grid.Ny-(pmlY+2):grid.Ny) = 1; %R
    %Q(grid.Nx-(pmlX+2):grid.Nx,:) = 1; %Z
         
    Q = diag(sparse(Q(:)));
    
    fsrc  = fsrc(:);
        
    A = Sx*DhZ/urr*Sx*DeZ;% + Sy*DhR/uzz*Sy*DeR + Sy*DhR/uzz*invRho;
    Az = Sx*DhZ/urr*Sx*DeZ;
    Ar = Sy*DhR/uzz*Sy*DeR;
    src = (Q*Ar-Ar*Q)*fsrc;
    Psi = A\src;
    Psi = reshape(Psi,grid.Nx,grid.Ny);
      imagesc(real(Psi));
      colorbar;
    