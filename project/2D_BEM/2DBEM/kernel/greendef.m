function [G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k,pxyb,betaP,xq,yq);

% [G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k,pxyb,betaP,xq,yq);
%
% Calculates the 2D Green's function for free field, rigid plane or
% plane with impedance.
%
% Note:
%   To obtain the incident pressure on the nodes from a line source in the absence
%   of the body using this function, do like this:
%
%    [G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k,position,betaP,xyb(:,1),xyb(:,2));
%    inc_pressure=i/4*(G0dir+G0ref)+Pbeta;
%
%   where 'position' are the [x,y] coordinates of the line source and xyb is the node
%   positions geometry file. 
%
%
% Input variables:
%   -k:       wavenumber.
%   -pxyb:    real vector containing the (x,y,body) values for
%             the point 'P', source or collocation point.
%   -betaP:    admitance of the plane. If it is NaN (not a number), the
%             free-field Green's function is calculated. If it is 0, the
%             plane is considered rigid.
%   -xq:      real column vector containing the global x-coordinates
%             for each point to calculate (e.g. integration points). 
%   -yq:      real column vector containing the global y-coordinates
%             for each point to calculate (e.g. integration points). 
%
%  Output variables:
%   -G0dir:    direct Green's function for free-field or rigid plane.
%   -G0ref:    reflected Green's function for free-field or rigid plane.
%   -dG0dirdR1:derivative of the direct Green's function.
%   -dG0dirdR2:derivative of the reflected Green's function.
%   -Pbeta:   correction term for a plane with finite impedance.
%   -dPbetadx,dPbetady :derivatives of the correction term.
%
% Reference:  -S.N. Chandler-Wilde and D.C. Hothersall: "Efficient calculation
%             of the Green function for acoustic propagation above an homogeneous
%             impedance plane", Journal of Sound and Vibration (1995) 180(5),
%             pp. 705-724.

% Vicente Cutanda 5-2001.


% Number of points on which to calculate the Green function (e.g., integration points).
npoints=length(xq);

% Gauss-Laguerre order
nlaguerre=15;
[t,w]=laguerre(nlaguerre,-0.5);  % basepoints(t) and weights(w) in Laguerre Quadrature

R1=sqrt((xq-pxyb(1)).^2+(yq-pxyb(2)).^2); % Distances from IPs to collocation/source point

R2=sqrt((xq-pxyb(1)).^2+(yq+pxyb(2)).^2); % Distances from IPs to image collocation/source point


% Free field Green function and its derivative
G0dir=besselh(0,1,k*R1);
dG0dirdR1=besselh(1,1,k*R1);


if isnan(betaP) % Free field case
   G0ref=zeros(npoints,1);
   dG0refdR2=zeros(npoints,1);

   Pbeta=zeros(npoints,1);
   dPbetadx=zeros(npoints,1);
   dPbetady=zeros(npoints,1);
   
elseif betaP==0 % Plane with infinite impedance
   G0ref=besselh(0,1,k*R2);
   dG0refdR2=besselh(1,1,k*R2);

	Pbeta=zeros(npoints,1);
   dPbetadx=zeros(npoints,1);
   dPbetady=zeros(npoints,1);
   
else % Plane with finite impedance
   G0ref=besselh(0,1,k*R2);
   dG0refdR2=besselh(1,1,k*R2);

   % Calculation of Pbeta and its derivative
   ro=k*R2.';
   gamma=repmat(((yq+pxyb(2))./R2).',nlaguerre,1);

   points=repmat(t,1,length(ro))./repmat(ro,nlaguerre,1);

   f=-(betaP+gamma.*(1+i*points))./(sqrt(points-i*2).*(points.^2-i*2*(1+betaP*gamma).*points-(betaP+gamma).^2));
	a=1+betaP*gamma-sqrt(1-betaP^2)*sqrt(1-gamma.^2);

   if betaP==1
     	cons=betaP*exp(i*ro)./(pi*sqrt(ro));
     	pgamma=-(cons.*w.')*f;
   
      if (imag(betaP)<0)&(real(a)>0)
        	ps=-betaP*exp(-i*ro*(1-a))/sqrt(1-betaP^2);
      elseif (imag(betaP)<0)&(real(a)==0)
        	ps=-betaP*exp(-i*ro*(1-a))/(2*sqrt(1-betaP^2));
      else
        	ps=zeros(size(ro));
      end
      	Pbeta=pgamma.'+ps.';
   else
      cons=betaP*exp(i*ro)./(pi*sqrt(ro));
    	g=f-(exp(-i*pi/4))*sqrt(a)./(2*sqrt(1-betaP^2).*(points-i*a));
     	P1=cons.*(w.'*g);
     	argument=exp(-i*pi/4)*sqrt(ro).*sqrt(a(1,:));
     	ercomplex=exp((argument/-i).^2).*eerfcc(argument/-i);
  	   P2=(betaP*exp(i*ro.*(1-a(1,:)))/(2*sqrt(1-betaP^2))).*ercomplex;
      Pbeta=-P1.'-P2.';   
   end
   
   if betaP==1
    	cons=(i*k*betaP*exp(i*ro)./(pi*sqrt(ro))).*sqrt(1-gamma(1,:).^2).*sign(xq.'-pxyb(1));
   
      f2=-(gamma+betaP.*(1+i*points))./(sqrt(points-i*2).*(points.^2-i*2*(1+betaP*gamma).*points-(betaP+gamma).^2));
  	   dpgamma=(cons.*w.')*f2;   
     
  	   if (imag(betaP)<0)&(real(a)<0)
     	   dps=i*k*betaP*sign(xq.'-pxyb(1)).*exp(i*ro.*(1-a(1,:)));
     	elseif (imag(betaP)<0)&(real(a)==0)
        	dps=0.5*i*k*betaP*sign(xq.'-pxyb(1)).*exp(i*ro.*(1-a(1,:)));
  	   else
         dps=zeros(size(ro));;
     	end
         dPbetadx=-dpgamma+dps;
   else
   
     	f2=-(gamma+betaP*(1+i*points))./(sqrt(points-i*2).*(points.^2-i*2*(1+betaP*gamma).*points-(betaP+gamma).^2));
     	cons=-i*k*betaP*exp(i*ro).*sign(xq.'-pxyb(1));
     	g=sqrt(1-gamma.^2).*f2-exp(-i*pi/4)*sqrt(a)./(2*(points-i*a));

      int=(1./(pi*sqrt(ro))).*(w.'*g);
      int2=0.5*exp(-i*ro.*a(1,:)).*ercomplex;

      dPbetadx=cons.*(int+int2);
   end

   dPbetady=(k*betaP/2)*besselh(0,1,ro)-i*k*betaP*Pbeta.';
   dPbetadx=dPbetadx.';
   dPbetady=dPbetady.';

end
