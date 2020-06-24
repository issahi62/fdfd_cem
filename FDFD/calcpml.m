function [ sx,sy ] = calcpml(grid,pmlX,pmlY)
    
    %Constants
    a_max = 3;
    sig_max = 1;
    p = 3;
    eta = 376.73032165;
    
    sx = ones(grid.Nx,grid.Ny);
    sy = ones(grid.Nx,grid.Ny);
    
    for n = 1 : pmlX
        % compute PML value
        a = 1 + a_max*(n/pmlX)^p;
        c = sig_max*sin(0.5*pi*n/pmlX)^2;
        s  = a*(1-1i*c*eta/(2*pi/grid.lam0));
        % incorporate value into PML      
        sx(pmlX-n+1,:) = s;
        sx(grid.Nx - pmlX + n,:) = s;
    end
    
    for n = 1 : pmlY
        % compute PML value
        a = 1 + a_max*(n/pmlY)^p;
        c = sig_max*sin(0.5*pi*n/pmlY)^2;
        s  = a*(1-1i*c*eta/(2*pi/grid.lam0));
        % incorporate value into PML      
        sy(:,pmlY-n+1) = s;
        sy(:,grid.Ny - pmlY + n) = s;
    end
end

