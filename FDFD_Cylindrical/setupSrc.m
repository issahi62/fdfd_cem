function [src] = setupSrc(grid,A,pmlX,pmlY,srcAngle)
degrees = pi / 180;
theta = srcAngle * degrees;
kinc = 2*pi / grid.lam0 * [sin(theta);cos(theta)];

% COMPUTE SOURCE FIELD
xa    = [0:grid.Nx-1]*grid.dx;
ya    = [0:grid.Ny-1]*grid.dy;
[Y,X] = meshgrid(ya,xa);

%fsrc  = exp((-1i*(kinc(1)*X))).*besselj(0,kinc(2)*Y);
fsrc = exp((-1i*(kinc(1)*X + kinc(2)*Y)));

fsrc  = fsrc(:);

% COMPUTE SCATTERED-FIELD MASKING MATRIX
Q = zeros(grid.Nx,grid.Ny);

Q(:,1:pmlY+2) = 1;
Q(1:pmlX+2,:) = 1;
Q(:,grid.Ny-(pmlY+2):grid.Ny) = 1;
Q(grid.Nx-(pmlX+2):grid.Nx,:) = 1;

Q = diag(sparse(Q(:)));

%Just TF useful debugging tool
% fsrcTF = (eye(grid.Nx*grid.Ny) - Q)*fsrc;
 %src = fsrcTF;
%   fsrcTF = reshape(fsrcTF,[grid.Nx,grid.Ny]);
%   imagesc(real(fsrcTF)')
%   colorbar;
%  fsrcSF = Q*fsrc;
%  fsrcSF = reshape(fsrcSF,[grid.Nx,grid.Ny]);
%imagesc(real(fsrcSF + fsrcTF)')

% COMPUTE SOURCE VECTOR
src = (Q*A-A*Q)*fsrc;
%src = fsrc;
%grid.Nx
%grid.Ny
%srcPlot = real(reshape(src,[grid.Nx,grid.Ny]))';
%imagesc(srcPlot)
%imagesc(real(src)')

end