 % Finite difference equations for cylinder and sphere
    % for 1D transient heat conduction with convection at surface
    % general equation is:
    % 1/alpha*dT/dt = d^2T/dr^2 + p/r*dT/dr for r ~= 0
    % 1/alpha*dT/dt = (1 + p)*d^2T/dr^2     for r = 0
    % where p is shape factor, p = 1 for cylinder, p = 2 for sphere
    function T = funcACbar(pbar,cpbar,kbar,h,Tinf,b,m,dr,dt,T)
    alpha = kbar./(pbar.*cpbar);    % effective thermal diffusivity
    Fo = alpha.*dt./(dr^2);         % effective Fourier number
    Bi = h.*dr./kbar;               % effective Biot number
    % [A] is coefficient matrix at time level n+1
    % {C} is column vector at time level n
    A(1,1) = 1 + 2*(1+b)*Fo(1);
    A(1,2) = -2*(1+b)*Fo(2);
    C(1,1) = T(1);
    for k = 2:m-1
        A(k,k-1) = -Fo(k-1)*(1 - b/(2*(k-1)));   % Tm-1
        A(k,k) = 1 + 2*Fo(k);                    % Tm
        A(k,k+1) = -Fo(k+1)*(1 + b/(2*(k-1)));   % Tm+1
        C(k,1) = T(k);
    end
    A(m,m-1) = -2*Fo(m-1);
    A(m,m) = 1 + 2*Fo(m)*(1 + Bi(m) + (b/(2*m))*Bi(m));
    C(m,1) = T(m) + 2*Fo(m)*Bi(m)*(1 + b/(2*m))*Tinf;
    % solve system of equations [A]{T} = {C} where temperature T = [A]\{C}
    T = A\C;
    end
