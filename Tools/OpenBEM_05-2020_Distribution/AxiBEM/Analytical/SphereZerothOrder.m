function p = SphereZerothOrder(k,c,rho,a,up,rtp);

% p = SphereZerothOrder(k,c,rho,a,up,rtp);
%
% Pulsating sphere vibrating on its zeroth order mode (U0)
% See for example Lecture Note "Radiation of Sound", by  Finn Jacobsen 
% and Peter M. Juhl. See also Pierce.
%
% Input:
%   -k:        Wavenumber, m-1. Can be a vector.
%   -c:        Speed of sound, m/s
%   -a:        Sphere radius, m
%   -rho:      Density of air, kg/m3
%   -up:       Maximum normal velocity on the sphere, m/s
%   -rtp:      Points where the sound pressure is calculated. They are 
%              external to the sphere. It is a column vector of r values.
%              There is no other dependance in this case.
%
% Output:
%   -pf:       Calculated sound pressure, Pa. Matrix with one row per 
%              point and one column per wavenumber.
%

% Vicente Cutanda Henriquez 05-2007

if any(rtp<a)
   error('There are points defined inside the sphere')
end

tmp=zeros(1,length(k));tmp(1,:)=k(1:end);k=tmp;clear tmp % k becomes a row vector, in any case
tmp=zeros(length(rtp),1);tmp(:,1)=rtp(1:end);rtp=tmp;clear tmp % rtp becomes a column vector, in any case

A0=ones(size(rtp))*(j*rho*c*a^2*up*k.*exp(j*k*a)./(1 + j*k*a));

p=A0.*exp(-j*rtp*k)./(rtp*ones(size(k))); % FJ (jwt convention)

%figure; plot(abs(p),'-x');grid
