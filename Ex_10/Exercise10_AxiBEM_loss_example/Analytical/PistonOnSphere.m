function p = PistonOnSphere(f,c,rho,a,alpha,up,R,theta);

% function p = PistonOnSphere(f,c,rho,alpha,up,R,theta);
%
% It Calculates the complex sound pressure (Pa) emitted by a cilindrcal
% vibrating piston on a spherical baffle.
%
%   f         vector of frecuencies (Hz)
%   c         sound speed (m/s)
%   rho       medium density  (kg/m³)
%   a         sphere radius (m)
%   alpha     piston angle (rad)
%   up        piston velocity amplitude (m/s)
%   R         measuring point radius (m)
%   theta     measuring point azimuth (rads)

% EARLY CALCULATIONS

k=2*pi.*f/c;            % wave number

% PERFORMING THE SUMMATORY
A=0;
for m=0:1:100
    A = A + P(m,theta) * (P(m-1,alpha)-P(m+1,alpha)) * (dhm(f,m,a,c)).^(-1).*hank(f,m,R,c);
end
% and calculating the sound pressure
p = - j*rho*c*up *A/2;

% SUBFUNCTIONS
% Legendre polynomial
function[P]=P(m,alpha)
if m==-1 P=1;
else P=legendre(m,cos(alpha)); end
P=P(1);

%Hankel function
function [hm]=hank(f,m,r,c) 

k=2*pi.*f/c;            % wave number
kr=k*r;
nu=m+(1/2);
hm=0;
H =0;
H = besselh(nu,2,kr);
hm=sqrt(pi./(2*kr)).*H;

%Derivative of Hankel function
function[dhm]=dhm(f,m,r,c) 
hm_1=0;
hm_2=0;
dhm=0;
hm_1=hank(f,m-1,r,c);
hm_2=hank(f,m+1,r,c);
dhm=(1/(2*m+1))*(m*hm_1-(m+1)*hm_2);
