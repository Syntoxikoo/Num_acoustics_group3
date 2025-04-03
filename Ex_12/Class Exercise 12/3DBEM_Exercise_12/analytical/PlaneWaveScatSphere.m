function [ptot,accuracy]=PlaneWaveScatSphere(k,a,r,theta,Ieps)

% ptot=PlaneWaveScatSphere(ka,r,theta,Ieps)
%
% Scattering by rigid sphere
% Bowman formulation
%
% Input:
%   -k:        Wavenumber, m-1
%   -a:        Sphere radius, m
%   -r:        Radii of the points where the sound pressure is calculated. 
%              They are external to the sphere. It is a vector.
%   -theta:    Theta angles of the points where the sound pressure is
%              calculated. It is a vector.
%   -Ieps:     Desired maximum error.
%
% Output:
%   -ptot:     Calculated sound pressure, Pa. Matrix with one row per 
%              r-coordinate and one column per theta-coordinate.
%

% Input:
% r=r/a radii normalized to sphere radius

tmp=zeros(1,length(theta));tmp(1,:)=theta(1:end);theta=tmp;clear tmp % theta becomes a row vector, in any case
tmp=zeros(length(r),1);tmp(:,1)=r(1:end);r=tmp;clear tmp % r becomes a column vector, in any case

ka=k*a;
kr=k*r;
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
