function [src] = setupSrc(grid,A,pmlX,pmlY)
degrees = pi / 180;
theta = 45 * degrees;
kinc = 2*pi / grid.lam0 * [sin(theta);cos(theta)];
%kinc = 1 / grid.lam0 * [sin(theta);cos(theta)];

% COMPUTE SOURCE FIELD
xa    = [0:grid.Nx-1]*grid.dx;
ya    = [0:grid.Ny-1]*grid.dy;
[Y,X] = meshgrid(ya,xa);
NumLambda = (2*pi);
fsrc  = exp(NumLambda / (2*pi) * (-1i*(kinc(1)*X+kinc(2)*Y)));

%fsrc = cos(kinc(1) * X) + cos(kinc(2) * Y);

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

% COMPUTE SOURCE VECTOR
src = (Q*A-A*Q)*fsrc;
%src = fsrcTF;
%grid.Nx
%grid.Ny
%srcPlot = real(reshape(src,[grid.Nx,grid.Ny]))';
%imagesc(srcPlot)
%imagesc(real(src)')

end