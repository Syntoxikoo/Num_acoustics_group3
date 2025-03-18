function p = SphereFirstOrder(k,c,rho,a,up,rtp);

% p = SphereFirstOrder(k,c,rho,a,up,rtp);
%
% Sphere vibrating on its first order mode (U0*cos(theta))
% See for example Lecture Note "Radiation of Sound", by  Finn Jacobsen 
% and Peter M. Juhl. See also Pierce.
%
% Input:
%   -k:        Wavenumber, m-1. Can be a vector.
%   -c:        Speed of sound, m/s
%   -a:        Sphere radius, m
%   -rho:      Density of air, kg/m3
%   -up:       Maximum normal velocity on the sphere, m/s
%   -rtp:      Points where the sound pressure is calcualted. They are 
%              external to the sphere. It is a matrix with two columns, 
%              r and theta coordinates of the points, m
%
% Ourput:
%   -pf:       Calculated sound pressure, Pa. Matrix with one row per 
%              point and one column per wavenumber.
%

% Vicente Cutanda Henriquez 05-2007

if any(rtp(:,1)<a)
   error('There are points defined inside the sphere')
end

tmp=zeros(1,length(k));tmp(1,:)=k(1:end);k=tmp;clear tmp % k becomes a row vector, in any case

A1=-ones(size(rtp,1),1)*(rho*c*a*up*k.*exp(j*k*a)./(1 - 2./(k*a).^2 + 2./(j*k*a))); 

p=(A1.*(1 + 1./(j*rtp(:,1)*k) )).*(exp(-j*rtp(:,1)*k)./(rtp(:,1)*k)).*(cos(rtp(:,2))*ones(1,length(k))); % FJ (jwt convention)


%figure; plot(abs(p),'-x');grid
