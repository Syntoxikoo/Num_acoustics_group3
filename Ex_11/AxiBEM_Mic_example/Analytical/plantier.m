function [pa,ph,var,vaz,vhr,vhz,vtr,vtz,er]=plantier(f,Rad,Hgh,r,z,n,pinc,Tm,mu_m,varargin);

% [pa,ph,var,vaz,vhr,vhz,vtr,vtz,er]=plantier(f,Rad,Hgh,r,z,n,pinc,Tm,mu_m,fact);
%
% Calculates the analytical solution for the sound field inside a thin cylindrical cavity
% with a coupled membrane and visco-thermal losses, as described in the reference.
%
% Input variables:
%   -f:      frequency, vector with frequencies to calculate.
%   -Rad:    radius of the cavity.
%   -Hgh:    heigth of the cavity.
%   -r:      values of the r-coordinate to calculate. Row vector.
%   -z:      values of the z-coordinate to calculate. Row vector.
%   -n:      number of terms in the expansions.
%   -pinc:   uniform incident pressure on the exterior of the membrane.
%   -Tm:     tension of the membrane.
%   -mu_m:   mass per unit area of the membrane.
%   -fact:   reduction factor of the viscisity (optional).
%            If fact > 1, the viscosity constants are reduced.
%            Suggested value fact = 100, to get rid of most viscosity effects.
%
% Output variables. Arranged in 3D matrices (except 'er'), the dimensions are respectively
% r-values, z-values and frequency values:
%   -pa:     acoustic pressure.
%   -ph:     thermal pressure.
%   -var:    radial component of the acoustic velocity.
%   -vaz:    normal component of the acoustic velocity.
%   -vhr:    radial component of the thermal velocity.
%   -vhz:    normal component of the thermal velocity.
%   -vtr:    radial component of the viscous velocity.
%   -vtz:    normal component of the viscous velocity.
%   -er:     normal displacement of the membrane, vector, one column per 'r', one row per frequency.
%
% Reference:
%   G. Plantier and M. Bruneau, "Heat conduction effects on the acoustic response of a membrane
%   separated by a very thin air film from a backing electrode". J. Acoustique 3 (1990), pp 243-250.

% Vicente Cutanda 5-2001

if nargin > 9
   fact=-abs(varargin{1});
else
   fact=-1;
end

% visco-thermal constants and convert to Plantier's notation.
[rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h]=VTconst(f,fact);
Lambda_a=1./tau_a;Lambda_h=1./tau_h;
alfa_a=phi_a./tau_a;alfa_h=phi_h./tau_h;
w=2*pi*f;K2=mu_m*w.^2/Tm;

% calculate zeros of Bessel function order 0.
xx=2;
kn=fzero('besselj(0,x)',xx)/Rad;
kni=kn;
for nn=2:n
   while round(kni*1e3)==round(kn(end)*1e3)
      xx=xx+1;
      kni=fzero('besselj(0,x)',xx)/Rad;
   end
   kn=[kn;kni];
end

kazn=sqrt(repmat(ka.^2,n,1)-repmat(kn.^2,1,length(f)));
khzn=sqrt(repmat(kh.^2,n,1)-repmat(kn.^2,1,length(f)));
kvzn=sqrt(repmat(kv.^2,n,1)-repmat(kn.^2,1,length(f)));

% calculate constants for all n-terms and frequencies and make the summations
pa=zeros(length(r),length(z),length(f));
ph=pa;var=pa;vaz=pa;vhr=pa;vhz=pa;vtr=pa;vtz=pa;
er=pinc./(Tm*repmat(K2.',1,length(r))).*(1-besselj(0,sqrt(K2).'*r)./besselj(0,repmat((sqrt(K2)*Rad).',1,length(r))));

zeros(length(f),length(r));
for ff=1:length(f)
   AB=zeros(6);RHS=zeros(6,1);
   for nn=1:n
      % build the system of equations and solve it to get a set of constants
      AB(1,[1 3])=1;
      AB(2,1:4)=[cos(kazn(nn,ff)*Hgh) sin(kazn(nn,ff)*Hgh) cos(khzn(nn,ff)*Hgh) sin(khzn(nn,ff)*Hgh)];
      AB(3,[2 4 5])=[alfa_a(ff)*kazn(nn,ff) alfa_h(ff)*khzn(nn,ff) 1];
      AB(4,[1 3 6])=[alfa_a(ff)*kn(nn) alfa_h(ff)*kn(nn) kvzn(nn,ff)/kn(nn)];
      AB(5,:)=[-alfa_a(ff)*kn(nn)*cos(kazn(nn,ff)*Hgh) -alfa_a(ff)*kn(nn)*sin(kazn(nn,ff)*Hgh) ...
               -alfa_h(ff)*kn(nn)*cos(khzn(nn,ff)*Hgh) -alfa_h(ff)*kn(nn)*sin(khzn(nn,ff)*Hgh) ...
               kvzn(nn,ff)/kn(nn)*sin(kvzn(nn,ff)*Hgh) -kvzn(nn,ff)/kn(nn)*cos(kvzn(nn,ff)*Hgh)];
      AB(6,:)=[(i*w(ff)*Lambda_a(ff)*cos(kazn(nn,ff)*Hgh))./(Tm*(kn(nn)^2-K2(ff)))+alfa_a(ff)*kazn(nn,ff)*sin(kazn(nn,ff)*Hgh) ...
               (i*w(ff)*Lambda_a(ff)*sin(kazn(nn,ff)*Hgh))./(Tm*(kn(nn)^2-K2(ff)))-alfa_a(ff)*kazn(nn,ff)*cos(kazn(nn,ff)*Hgh) ...
               (i*w(ff)*Lambda_h(ff)*cos(khzn(nn,ff)*Hgh))./(Tm*(kn(nn)^2-K2(ff)))+alfa_h(ff)*khzn(nn,ff)*sin(khzn(nn,ff)*Hgh) ...
               (i*w(ff)*Lambda_h(ff)*sin(khzn(nn,ff)*Hgh))./(Tm*(kn(nn)^2-K2(ff)))-alfa_h(ff)*khzn(nn,ff)*cos(khzn(nn,ff)*Hgh) ...
               -cos(kvzn(nn,ff)*Hgh) -sin(kvzn(nn,ff)*Hgh)];
      RHS(6)=2*i*w(ff)*pinc/(Tm*(kn(nn)^2-K2(ff))*kn(nn)*Rad*besselj(1,kn(nn)*Rad));
      const=AB\RHS;
      Aan=const(1);Ban=const(2);Ahn=const(3);Bhn=const(4);Avn=const(5);Bvn=const(6);
      % calculate acoustic variables
      pa(:,:,ff)=pa(:,:,ff) + besselj(0,kn(nn)*r).'*(Aan*cos(kazn(nn,ff)*z)+Ban*sin(kazn(nn,ff)*z))*Lambda_a(ff);
      ph(:,:,ff)=ph(:,:,ff) + besselj(0,kn(nn)*r).'*(Ahn*cos(khzn(nn,ff)*z)+Bhn*sin(khzn(nn,ff)*z))*Lambda_h(ff);
      var(:,:,ff)=var(:,:,ff) + kn(nn)*besselj(1,kn(nn)*r).'*(Aan*cos(kazn(nn,ff)*z)+Ban*sin(kazn(nn,ff)*z))*alfa_a(ff);
      vaz(:,:,ff)=vaz(:,:,ff) - besselj(0,kn(nn)*r).'*(Aan*sin(kazn(nn,ff)*z)-Ban*cos(kazn(nn,ff)*z))*alfa_a(ff)*kazn(nn,ff);
      vhr(:,:,ff)=vhr(:,:,ff) + kn(nn)*besselj(1,kn(nn)*r).'*(Ahn*cos(khzn(nn,ff)*z)+Bhn*sin(khzn(nn,ff)*z))*alfa_h(ff);
      vhz(:,:,ff)=vhz(:,:,ff) - besselj(0,kn(nn)*r).'*(Ahn*sin(khzn(nn,ff)*z)-Bhn*cos(khzn(nn,ff)*z))*alfa_h(ff)*khzn(nn,ff);
      vtr(:,:,ff)=vtr(:,:,ff) + kvzn(nn,ff)/kn(nn)*besselj(1,kn(nn)*r).'*(Avn*sin(kvzn(nn,ff)*z)-Bvn*cos(kvzn(nn,ff)*z));
      vtz(:,:,ff)=vtz(:,:,ff) + besselj(0,kn(nn)*r).'*(Avn*cos(kvzn(nn,ff)*z)+Bvn*sin(kvzn(nn,ff)*z));
      er(ff,:)=er(ff,:) + 1/(Tm*(kn(nn)^2-K2(ff)))*besselj(0,kn(nn)*r)*...
                          (Lambda_a(ff)*(Aan*cos(kazn(nn,ff)*Hgh)+Ban*sin(kazn(nn,ff)*Hgh))+...
                          Lambda_h(ff)*(Ahn*cos(khzn(nn,ff)*Hgh)+Bhn*sin(khzn(nn,ff)*Hgh)));
   end
end
