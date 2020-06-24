function [src] = setupSrc(grid,A,pmlX,pmlY,srcAngle)
degrees = pi / 180;
theta = srcAngle * degrees;
src.kinc = 2*pi / grid.lam0 * [sin(theta);cos(theta)];
%kinc = 1 / grid.lam0 * [sin(theta);cos(theta)];

% COMPUTE SOURCE FIELD
xa    = [0:grid.Nx-1]*grid.dx;
ya    = [0:grid.Ny-1]*grid.dy;
[Y,X] = meshgrid(ya,xa);
NumLambda = 1;
%fsrc  = exp(NumLambda * (-1i*(src.kinc(1)*X+src.kinc(2)*Y)));
w = 8;
k0 = 2*pi / grid.lam0;
fsrc = exp(-(k0*(X-2.5*grid.lam0)/w).^2) .* exp(1j * k0 *Y); 
src.fsrc = fsrc;
subplot(2,2,2)
imagesc(real(fsrc)');
colorbar;
title('Source');

fsrc  = fsrc(:);

% COMPUTE SCATTERED-FIELD MASKING MATRIX
Q = zeros(grid.Nx,grid.Ny);

Q(:,1:pmlY+2) = 1;
Q(1:pmlX+2,:) = 1;
Q(:,grid.Ny-(pmlY+2):grid.Ny) = 1;
Q(grid.Nx-(pmlX+2):grid.Nx,:) = 1;

% Q(:,1:20+2) = 0;
% Q(1:20+2,:) = 0;
% Q(:,grid.Ny+1-(20+2):grid.Ny) = 0;
% Q(grid.Nx+1-(20+2):grid.Nx,:) = 0;

Q = diag(sparse(Q(:)));

%Just TF useful debugging tool
% fsrcTF = (eye(grid.Nx*grid.Ny) - Q)*fsrc;
 %src = fsrcTF;
%  fsrcTF = reshape(fsrcTF,[grid.Nx,grid.Ny]);
%  imagesc(real(fsrcTF)')
%  fsrcSF = Q*fsrc;
%  fsrcSF = reshape(fsrcSF,[grid.Nx,grid.Ny]);
%imagesc(real(fsrcSF + fsrcTF)')

src.SF = reshape(Q*fsrc,[grid.Nx,grid.Ny]);
src.SF(1:pmlX+1,:) = 0;

% COMPUTE SOURCE VECTOR
src.vec = (Q*A-A*Q)*fsrc;
%src = fsrcTF;
%grid.Nx
%grid.Ny
%srcPlot = real(reshape(src,[grid.Nx,grid.Ny]))';
%imagesc(srcPlot)
%imagesc(real(src)')

end