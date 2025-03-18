function pana=ana_sphere(ka,xyzb,Ieps);
%analyt=ana_sphere(ka,xyzb,Ieps);
%Computes the analytical solution for scattering by a rigid sphere
%of a plane wave
%
%ka - dimensionless frequency (wavenumber times radius)
%xyzb - each row in xyzb contains calculation point
%Iteration epsilon

theta=acos(xyzb(:,3));
pana= sphere(ka,1,theta,Ieps)';


function [ptot,accuracy]=sphere(ka,r,theta,Ieps)
% ptot=sphere(ka,r,theta,Ieps)
% r=r/a radii normalized to sphere radius
% r must be a column
% theta must be a row
% scattering by rigid sphere
% Bowman formulation

kr=ka.*r;
n=0;
ptot=zeros(length(r),length(theta));
accuracy=[];
rel_accuracy=Ieps*10;
while rel_accuracy>Ieps
   apri = (n*besselj(n-0.5,ka)-(n+1)*besselj(n+1.5,ka))./(n*besselh(n-0.5,ka)-(n+1)*besselh(n+1.5,ka));
   Pmatrix=legendre(n,cos(theta));
   P=Pmatrix(1,:); % vector of legendre_n function of different arguments
   factor=(-i)^n*(2*n+1)*sqrt(pi/2./kr);
   jn=besselj(n+0.5,kr);
   hn=besselh(n+0.5,kr);
   iter=(factor.*(jn-apri.*hn))*P;
   ptot=ptot+iter;
   [rel_accuracy]=min(max(abs(iter./ptot)'));
   [accuracy]=[accuracy; rel_accuracy];
   n=n+1;
end


