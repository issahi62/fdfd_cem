clear;
clc;
sensor_location = [102;102];
sensor_location = 6.8359e-5 * sensor_location;
rho   = norm(sensor_location);
krho = 8.716;
phi = angle(sensor_location(1)+1j*sensor_location(2));
 nu    = -50:50; % order of the Bessel function
 x   = 2.734;
 a_n = -besselj(nu,x)./besselh(nu,2,x);
 E_z = sum(1j.^(-nu).*a_n.*besselh(nu,2,krho).*exp(1j*nu*phi));
 E_z