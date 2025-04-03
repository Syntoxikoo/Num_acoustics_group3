function [p, v_r, v_theta, v_rA, v_thetaA, v_rV, v_thetaV] = SphereFirstOrder(k,c,rho,a,up,rtp,varargin);

% [p,v_r,v_theta,v_rA,v_thetaA,v_rV,v_thetaV] = SphereFirstOrder(k,c,rho,a,up,rtp,S,kv);
%
% Sphere vibrating on its first order mode (up*cos(theta))
% See for example Lecture Note "Radiation of Sound", by  Finn Jacobsen 
% and Peter M. Juhl. See also Pierce.
% If kv is included, the solution has viscous losses. The reference is 
% section 6.9 in 'Elements of Acoustics' by Temkin.
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
%   -S:        sign convention (1, Karra exp(-jwt) or -1, Bruneau exp(jwt)).
%              If not supplied, -1 is used (exp(jwt) convention). 
%   -kv:       if this parameter exists, the solution includes viscous 
%              losses with a viscous wavenumber kv. Can be a vector,
%              following k.
%
% Ourput:
%   -p:        Calculated sound pressure, Pa. Matrix with one row per 
%              point and one column per wavenumber.
%   -v_r:      Total velocity in radial direction. Same size as p.
%   -v_theta:  Total velocity in theta direction. Same size as p.
%   -v_rA:     Acoustic (irrotational) velocity in radial direction. Same size as p.
%   -v_thetaA: Acoustic (irrotational) velocity in theta direction. Same size as p.
%   -v_rV:     Viscous (rotational) velocity in radial direction. Same size as p.
%   -v_thetaV: Viscous (rotational) velocity in theta direction. Same size as p.
%
%  The last four outputs are empty if kv is not supplied (lossless calculation)

% Peter Møller Juhl & Vicente Cutanda Henriquez 03-2012

if any(rtp(:,1)<a)
   error('There are points defined inside the sphere')
end
tmp=zeros(1,length(k));tmp(1,:)=k(1:end);k=tmp;clear tmp % k becomes a row vector, in any case

if nargin>6  % Sign convention 
    S=sign(varargin{1});
else
    S=-1;
end

v_rA=[];v_thetaA=[];v_rV=[];v_thetaV=[];

if nargin>7
    % calculation including losses (uses exp(-jwt) convention)
    kv=varargin{2};tmp=zeros(1,length(kv));tmp(1,:)=kv(1:end);kv=tmp;clear tmp % kv becomes a row vector, in any case
    kv=real(kv)+1j*abs(imag(kv)); % impose exp(-jwt) convention to the viscous wavenumber
    
    Ka=kv*a; % kv in VTconst correspond to K in Temkin
    ka=k*a; % use perfect fluid wavenumber due to remarks above (6.9.1)
    % Note that here ka is wavenumber times radius - not acoustic wavenumber
    % as obtained from VTconst
    
    % Setup and solve (6.9.14-15)
    [h1_ka, dh1_ka] = sph_hankel(ka);
    [h1_Ka, dh1_Ka] = exp_sph_hankel(Ka,Ka); % scale to avoid very small coef's
    a11=3*1j*k.*dh1_ka; a12=-6*1j/a*h1_Ka;
    a21=3*1j/a*h1_ka; a22=-3*1j/a*(Ka.*dh1_Ka+h1_Ka);
    
    A1=zeros(size(k)); exp_beta_B1=zeros(size(k));
    for kk=1:length(k)
        [AB]=[a11(kk) a12(kk); a21(kk) a22(kk)]\[up; up];
        A1(1,kk)=AB(1);
        exp_beta_B1(1,kk)=AB(2);
    end
    A1=ones(size(rtp,1),1)*A1;
    exp_beta_B1=ones(size(rtp,1),1)*exp_beta_B1;
    
    % Find pressure potential (6.9.2)
    kr=rtp(:,1)*k;
    [h1_kr,dh1_kr]=sph_hankel(kr);
    phi=1j*3*A1.*h1_kr.*(cos(rtp(:,2))*ones(1,length(k))); % should become a matrix
    p=1j*rho*c*(ones(size(rtp,1),1)*k).*phi; % (6.7.17)
    
    % Find radial velocity (6.9.11) with n=1
    Kr=rtp(:,1)*kv;
    [exp_h1_Kr,exp_dh1_Kr]=exp_sph_hankel(Kr,ones(size(rtp,1),1)*Ka); %remember the scale on B1
    v_r=1j*3*(ones(size(rtp,1),1)*k).*(A1.*dh1_kr-2*exp_beta_B1./kr.*exp_h1_Kr).*(cos(rtp(:,2))*ones(1,length(k)));
    
    v_rA=1j*3*(ones(size(rtp,1),1)*k).*(A1.*dh1_kr).*(cos(rtp(:,2))*ones(1,length(k)));
    v_rV=1j*3*(ones(size(rtp,1),1)*k).*(-2*exp_beta_B1./kr.*exp_h1_Kr).*(cos(rtp(:,2))*ones(1,length(k)));
    
    % Find theta component of velocity (6.9.13) with n=1
    v_theta=1j*3./(rtp(:,1)*ones(1,length(k))).*(A1.*h1_kr-exp_beta_B1.*(Kr.*exp_dh1_Kr+exp_h1_Kr)).*(-sin(rtp(:,2))*ones(1,length(k)));
    
    v_thetaA=1j*3./(rtp(:,1)*ones(1,length(k))).*(A1.*h1_kr).*(sin(rtp(:,2))*ones(1,length(k)));
    v_thetaV=1j*3./(rtp(:,1)*ones(1,length(k))).*(-exp_beta_B1.*(Kr.*exp_dh1_Kr+exp_h1_Kr)).*(sin(rtp(:,2))*ones(1,length(k)));
    
    
    if S==-1, % Change to exp(jwt) convention if required
        p=conj(p);
        v_r=conj(v_r); v_theta=conj(v_theta);
        v_rA=conj(v_rA); v_thetaA=conj(v_thetaA);
        v_rV=conj(v_rV); v_thetaV=conj(v_thetaV);
    end
    
else
    % calculation without losses  (uses exp(jwt) convention)
    A1=-ones(size(rtp,1),1)*(rho*c*a*up*k.*exp(1i*k*a)./(1 - 2./(k*a).^2 + 2./(j*k*a)));  % (5.38)
    
    p=(A1.*(1 + 1./(1i*rtp(:,1)*k) )).*(exp(-1i*rtp(:,1)*k)./(rtp(:,1)*k)).*(cos(rtp(:,2))*ones(1,length(k))); % FJ (jwt convention)
    v_r=-A1/(rho*c).*(1-2./(rtp(:,1)*k).^2+2./(1j*rtp(:,1)*k)).*(exp(-1j.*(rtp(:,1)*k))./(rtp(:,1)*k)).*(cos(rtp(:,2))*ones(1,length(k))); %(5.35)
    v_theta=(1j.*A1)/(rho*c).*(1+1./(1j*rtp(:,1)*k)).*(exp(-1j*rtp(:,1)*k)./((rtp(:,1).^2)*k)).*(sin(rtp(:,2))*ones(1,length(k))); % (5.36)
    
    if S==1, % Change to exp(-jwt) convention if required
        p=conj(p);
        v_r=conj(v_r); v_theta=conj(v_theta);
    end
    
end


end



function [h1,dh1]=sph_hankel(z)
% returns spherical hankel function 1st kind 1st order and its derivative
h1=exp(1j*z)./z.*(1./(1j*z)-1);
dh1=exp(1j*z)./z.*(-1j+2./z+2*1j./(z.^2));
end

function [exp_h1,exp_dh1]=exp_sph_hankel(z,zn)
% returns spherical hankel function 1st kind 1st order and its derivative
% multiplied with exp(-j*zn) in order to scale for large imaginary arguments
exp_h1=exp(1j*(z-zn))./z.*(1./(1j*z)-1);
exp_dh1=exp(1j*(z-zn))./z.*(-1j+2./z+2*1j./(z.^2));
end
