function [rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu]=VTconst(f,varargin)

% [rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu]=VTconst(f,S)
%
% *****************************
% *  VISCO-THERMAL CONSTANTS  *
% *****************************
%
% Calculates visco-termal constants as defined in the paper:
% 
% Bruneau M., Herzog Ph., Kergomard J and Polack J.D.
% "General formulation of the dispersion equation in bounded visco-thermal
% fluid, and application to some simple geometries"
% Wave Motion 11 (1989), pp. 441-451
%
% Input:
%    -f: vector with frequencies to be used.
%    -S: sign convention (1, Karra exp(-jwt) or -1, Bruneau exp(jwt)).
%        If not supplied, -1 is used (exp(jwt) convention). 
%        If abs(S) is not 1, it is a constant that divides the viscosity.
%        Suggested value S = -100, to get rid of most viscosity effects.
%
% Output: (all vectors of the same length as 'f')
%   -kp: perfect fluid wavenumber
%   -ka: acoustic wavenumber
%   -kh: thermal wavenumber
%   -kv: viscous wavenumber
%   -tau_a, tau_h, phi_a, phi_h: constants in Bruneau's model.
%   -eta:    2nd viscosity of air, N*s/m2 = Pa*s
%   -mu:     Viscosity of air, N*s/m2 = Pa*s


% Vicente Cutanda 11-2000

if nargin > 1
   S=sign(varargin{1});
   fact=abs(varargin{1});
else
   S=-1;
   fact=1;
end
% physical constants (MKS units) of air in normal conditions:
pa = 101325;          % Static pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Relative humidity (%)
[rho,c,cf,cpcv,mu,alfa,Cp,Cv,landa,beta]=amb2prop(pa,t,Hr,1000); % K. Rasmussen calls the dynamic viscosity coefficient "eta", here it is renamed to "mu" (VCH 6-2011)
%eta=110e-7;           % 2nd (or bulk) viscosity (N*s/m2)
eta=0.6*mu;           % 2nd (or bulk) viscosity (N*s/m2). Pierce eq. 10-7.23 
%[eta,mu,rho,c,cpcv,landa,Cp,Cv,beta]=PhysicalConst(fact); % old version

% characteristic lengths:
lh=landa/(rho*Cp*c);
lv=(eta+4*mu/3)/(rho*c);
lvp=mu/(rho*c);
lvh=lv+(cpcv-1)*lh;
lvhp=(cpcv-1)*(lh-lv);

% wavenumbers:
kp=2*pi*f/c;                              %perfect fluid
ka2=kp.^2./(1-S*i*kp*lvh-kp.^2*lh*lvhp);   %acoustic^2   Simplified: ka2=kp.^2.*(1+S*i*kp*lvh);
kh2=S*i*kp./(lh*(1+S*i*kp*lvhp));          %entropic^2   Simplified: kh2=S*i*kp.*(1-S*i*kp*lvhp)/lh;
kv2=S*i*kp/lvp;                           %viscous^2

% The result of the square root is considered positive (??):
ka=sqrt(ka2);   %acoustic
kh=sqrt(kh2);   %entropic
kv=sqrt(kv2);   %viscous

% thermal boundary condition constants:
tau_a=(cpcv-1)./(beta*cpcv*(1+S*i*lh*ka2./kp));
tau_h=(cpcv-1)./(beta*cpcv*(1+S*i*lh*kh2./kp));

% velocity boundary condition constants:
phi_a=-S*i./(rho*2*pi*f.*(1+S*i*lv*ka2./kp));
phi_h=-S*i./(rho*2*pi*f.*(1+S*i*lv*kh2./kp));

