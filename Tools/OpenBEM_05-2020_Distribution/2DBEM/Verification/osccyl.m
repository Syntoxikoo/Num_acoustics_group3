function p=osccyl(k,a,xyb)
% p=osccyl(k,a,xyb)
% computes the analytical solution for radiation of oscillating
% cylinder in 2-D
% k wavenumber (not yet vecorized)
% a radius

ka=k*a;
U1=1; % 1 m/s surface velocity at 0 degrees

% K. Rasmussen p. 173 for velocity and velocity potential
H12=besselh(1,2,ka);
H02=besselh(0,2,ka);
B1=a*U1/(H12-ka*H02);
phi=B1*xyb(:,1)*H12;

% change from exp(jwt) to exp(-jwt)
phi=conj(phi);

% the BEM code works with dp/dn and p
% for dp/dn=1 U1=-1/(jw rho) dp/dn
% which means that phi should be multiplied with -1/(jw rho)
% and to translate into the pressure p = jw*rho*phi
% all in all we find 
p=-phi;