    % Kinetic equations for reactions of wood, first-order, Arrhenious type     equations
    % K = A*exp(-E/RT) where A = pre-exponential factor, 1/s
    % and E = activation energy, kJ/mol
    
function [rww, rwg, rwt, rwc] = funcY(A1,E1,A2,E2,A3,E3,R,T,pww)
    K1 = A1.*exp(-E1./(R.*T));    % wood -> gas (1/s)
    K2 = A2.*exp(-E2./(R.*T));    % wood -> tar (1/s)
    K3 = A3.*exp(-E3./(R.*T));    % wood -> char (1/s)
    rww = -(K1+K2+K3).*pww;      % rate of wood consumption (rho/s)
    rwg = K1.*pww;               % rate of gas production from wood (rho/s)
    rwt = K2.*pww;               % rate of tar production from wood (rho/s)
    rwc = K3.*pww;               % rate of char production from wood (rho/s)
end