function [Bcontrib,Acontrib,CK]=intF2(przb,elknrzb,m,przn)

%  [F2A,F2B,CK]=intF2(przb,elknrzb,m)
%  
%  Calculates the singular part of the FA and FB integrals
%
%  Input:
%    przb:    real vector containing the (rho,z,body) values for
%             the point 'P'
%    elknrzb: real matrix, one row for each node in the element
%             each row contains (rho,z,body) for the node
%    m:       integer scalar - circumferential mode number
%    przn:    normal vector at the collocation point (viscous mode calculation)
%
%  Output:
%    F2A: complex vector, contains the contribution of each
%         node in the element to the B matrix
%
%    F2B: complex vector, contains the contribution of each
%         node in the element to the A matrix
%  
%    CK:  real vector, contains the contribution of each
%         node in the element to the C constants
%  

%  Revised pmj 990128 to include B matrix
%  Integration density may now be different for angular integration
%  By msj 990115
%
%  From 'integrnd.m' by Vicente Cutanda 1998
%  
% Based on P. Juhl
% 'An axisymmetric integral equation formulation for free space
%  non-axisymmetric radiation and scattering of a known incident wave'
% Journal of Sound and Vibration (1993) vol.163 pp.397-406

persistent COEFFS % These coefficients are constant in all elements
                  % They are quite expensive to evaluate
if isempty(COEFFS)
   for mm=0:24
      for ii=mm:-1:0
         n=2*(mm-ii);
         if mm==0
            c(n+1)=1;
         elseif n==0
            c(n+1)=((-1)^ii)*2;
         else
            c(n+1)=-c(n-1).*(mm+n/2-1).*(mm-n/2+1)/(n-1)/n;
         end
%         if ii~=0 % Debugging test
%            if c(n+1)~=((-1)^ii)*2*mm/ii*nchoosek(2*mm-ii-1,ii-1)
%              c(n+1) - ((-1)^ii)*2*mm/ii*nchoosek(2*mm-ii-1,ii-1), c(n+1), n
%              error('Error in COEFFS');
%            end
%         end
      end
      COEFFS(mm+1,1:2*mm+1)=c;
   end
end

[nknel,dummy]=size(elknrzb);

% Pre-check to see if the calculation point is close to the element, and
% obtain integration points accordingly
elknrz=elknrzb(:,1:2);prz(1,1:2)=przb(1:2);
maxnoddist=max(sqrt(sum(diff([elknrz; elknrz(1,:)]).^2,2))); % Get size of the maximum distance between nodes
maxdist=max(sqrt(sum((elknrz-ones(size(elknrz,1),1)*prz).^2,2))); % distance from the calculation point to the farthest node in the element:

if maxdist>maxnoddist*1.1 %& ~is_vt% it should not get nodes from adjoining elements
    [bp,wf]=rgauss(nknel-1,20);
else % interval division is used also for visco-thermal modes (VCH 7-2011)
    
    % The function nsingrule takes care of near-singular AND singular (diagonal
    % coefficients) integrands, by means of a limited recursive interval
    % division. It may be convenient to implement a specialized singular
    % numerical integration rule to solve the singular case.
%     if is_vt
% %        [bp,wf]=rgauss(nknel-1,20);
%         [bp,wf]=nsingrule(20,elknrzb,przb);
%     else
        [bp,wf]=nsingrule(8,elknrzb,przb);
%     end
end

[psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,bp);
jacobi=sqrt(nrho.^2+nz.^2);

x=bp;


% [nknel,dummy]=size(elknrzb);
% 
% %  In my opinion the integration order needed for the generator integral
% %  is determined by the order of shape function used when not dealing with (near)singular integrals (pmj)
% 
% [bp,wf]=gaussrule(10);
% x=bp;
% [psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
% % Note that 'nrho' and 'nz' have not been normalised, so
% % The Jacobian should only be used for 'F1A' and 'F2A' (msj)
% jacobi=sqrt(nrho.^2+nz.^2);
% 
% % Angular integration of singular part ('F2B' in thesis)
% % Done analytically with elliptic integrals
% 
% % First see if integration needs to be refined due to near singularity
% % or diagonal term
% threshold=0.2;
% dist2=(elknrzb(:,1)-przb(1)).^2+(elknrzb(:,2)-przb(2)).^2;
% ell2=(elknrzb(1,1)-elknrzb(nknel,1)).^2+(elknrzb(1,2)-elknrzb(nknel,2)).^2;
% [mindist2,closenode]=min(dist2);
% 
% if mindist2<ell2
%    % Point 'P' is closer than one element length from the element
%    % so it is necessary to check if the integration needs to be refined
%    dist2=(rhoq-przb(1)).^2+(zq-przb(2)).^2;
%    [mindist2q,closeq]=min(dist2);
% 
%    if mindist2q<mindist2
%       mindist2=mindist2q;
%       singpos=x(closeq);
%    else
%       singpos=closenode-2;
%    end
% 
%    if mindist2<ell2*threshold^2
%       if mindist2<ell2*10*eps^2
%          % Diagonal term
%          mindist2=0.0001;
%       else
%          % Near singularity
%          disp('Near Singular Integration Performed');
%       end
% 
%       [x,wf]=nsingrule(8,singpos,2*sqrt(mindist2/ell2));
% 
%       [psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
%       jacobi=sqrt(nrho.^2+nz.^2);
% 
%       if mindist2~=0.0001 & min((rhoq-przb(1)).^2+(zq-przb(2)).^2)<0.01*mindist2
%       % Yes it was a bit dirty of me to test for diagonal entries in this way
%          warning('Oblique collocation point too close to the boundary');
%          disp('PRESS RETURN TO CONTINUE');
%          pause;
%       end
%    end
% end

% Angular integration
dz2=(zq-przb(2)).^2;
rstr2=(rhoq+przb(1)).^2+dz2;
rstr=sqrt(rstr2);
%kstr2=4./rstr2.*rhoq*przb(1); % First estimation, can produce numerical errors
%kstr=sqrt(kstr2);

% The following formula avoids FP errors when used with Landen transformation (VC 6-2011)
kdash2=((rhoq-przb(1)).^2+dz2)./((rhoq+przb(1)).^2+dz2);
kdash=sqrt(kdash2);
kstr=sqrt(1-kdash2);   % VCH 6-2011
kstr2=(1-kdash2);   %  VCH 6-2011

if min(kstr2)>0.5  % Use elliptic integrals

   %  Elliptic integral (Matlab function, works on matrices):
   %[kea,ea]=ellipke(kstr2); % removed by VCH 6-2011

   %  Elliptic integral with descending Landen transformation, see Abramowitz or NIST handbook (VCH 6-2011):
   kstr_1=(1-kdash)./(1+kdash);
   [kea_1,ea_1]=ellipke(kstr_1.^2); % kstr_1 is computed more accurately than kstr  
   kea=(1+kstr_1).*kea_1; % complete elliptic integral of the first kind.
   ea=(1+kdash).*ea_1-kdash.*kea; % complete elliptic integral of the second kind.

   dkstrr=sqrt(przb(1)./rhoq)./rstr -kstr.*(rhoq+przb(1))./rstr2;
   dkstrz=-kstr.*(zq-przb(2))./rstr2;
   % eps^2 introduced to avoid division by zero
   dke=(ea-kdash2.*kea)./(kstr+eps^2)./kdash2;

   drstrr=(rhoq+przb(1))./rstr;
   drstrz=(zq-przb(2))./rstr;

   Xm=zeros(size(rstr));
   dXm=zeros(size(rstr));
   for ii=m:-1:0
      n=2*(m-ii);
      if abs(przb(1))<1e-6
         % Special case because kstr==0
         dIn=dke;
         if n==0
            In=kea;
         else
            Inm2=In;
            In=(n-1)/n.*Inm2;
         end
      else
         if n==0
            In=kea;
            dIn=dke;
         elseif n==2
            Inm2=In;
            In=(ea-kdash2.*kea)./kstr2;
            dInm2=dIn;
            dIn=((2-kstr2).*kea-2*ea)./(kstr.*kstr2);
         else
            Inm4=Inm2;
            Inm2=In;
            In=((n-2).*(2*kstr2-1).*Inm2+...
                (n-3).*kdash2.*Inm4)./((n-1)*kstr2);
            dInm4=dInm2;
            dInm2=dIn;
            dIn=((n-2).*((2*kstr2-1).*dInm2+2./kstr.*Inm2)+...
                (n-3).*(kdash2.*dInm4-2./kstr.*Inm4))./((n-1)*kstr2);
         end
      end
      if any(In<0) | any(In>pi/2./kdash)
         In, pi/2./kdash, n, m
         error('In<0 or In>pi/(2*kdash)');
      end
      I(:,n+1)=In; % I(1)=I0 since matlab arrays must start at 1
      dI(:,n+1)=dIn;
   end


   Xm=I(:,2*m+1);
   dXm=dI(:,2*m+1);
   for n=2*m-2:-2:0
      ii=m-n/2;
      Xm=Xm*4+((-1)^ii)*2*m/ii*nchoosek(2*m-ii-1,ii-1)*I(:,n+1);
      dXm=dXm*4+((-1)^ii)*2*m/ii*nchoosek(2*m-ii-1,ii-1)*dI(:,n+1);
   end
   if m~=0
      Xm=Xm./2;
      dXm=dXm./2;
   end
   
   F2B=(-1)^m*4./(rstr2).*((rstr.*dXm.*dkstrr-Xm.*drstrr).*nrho+(rstr.*dXm.*dkstrz-Xm.*drstrz).*nz);
   Cfunk1a=4./(rstr2).*((rstr.*dke.*dkstrr-kea.*drstrr).*nrho+(rstr.*dke.*dkstrz-kea.*drstrz).*nz);
   F2A=((-1)^m)*4./rstr.*Xm;

   % Correction of vector directions of P and Q  >>>> theta dependence not
   % yet implemented, to be done analytically?
   dotprodPQ=((nrho./jacobi)*przn(1)) + ((nz./jacobi)*przn(2)); 
   F2A=F2A.*dotprodPQ;
 %  F2B=F2B.*dotprodPQ;

else % Do the angular integration numerically
    
%     if is_vt
%         % When the wavenumber has a significant imaginary part (viscous and
%         % thermal modes), the integral is calculated using the near-singular
%         % basepoints and weights centered at theta=0. (VCH 7-2011)
%         [bp_ang,wf_ang]=nsingrule(8,[-1 0 1;0 0 1;1 0 1],[-1 0 1]);
%     else
        nip=8*m+8;
        bp_ang=[linspace(-1+1/nip,1-1/nip,nip)'];
        wf_ang=[2*ones(nip,1)./nip];
        if abs(sum(wf_ang)-2)>1e-4 | find(abs(bp_ang)>=1)
            sum(wf_ang), bp_ang
            error('Error in integration formula');
        end
%     end
    
    theta=bp_ang*(pi/2)+pi/2;
    bpones=ones(size(bp_ang));
    xones=ones(size(x));
    Ra=sqrt((rhoq.^2)*bpones'+przb(1)^2+((zq-przb(2)).^2)*bpones' ...
        -2*przb(1)*rhoq*cos(theta)');
    Bfr1=1./Ra;
    fr1=(-1./Ra.^3).* ...
        ((rhoq*bpones'-przb(1)*xones*cos(theta)').*(nrho*bpones')+ ...
        (zq-przb(2)).*nz*bpones');
    
    Cfunk1a=2*fr1*wf_ang*(pi/2);
    
    % Correction of vector directions of P and Q
    dotprodPQ=((nrho./jacobi)*przn(1))*cos(theta)' + ((nz./jacobi)*przn(2))*bpones'; 
%    fr1=fr1.*dotprodPQ;
    Bfr1=Bfr1.*dotprodPQ;
    
    F2B=2*fr1*(wf_ang.*cos(m*theta))*pi/2;
    F2A=2*Bfr1*(wf_ang.*cos(m*theta))*pi/2;
end

% % Correction of vector directions of P and Q  >>>> theta dependence not
% % yet implemented, to be done analytically?
dotprodPQ=((nrho./jacobi)*przn(1)) + ((nz./jacobi)*przn(2));
% F2A=F2A.*dotprodPQ;
F2B=F2B.*dotprodPQ;

% Integration over element length
Acontrib=[wf' *(((rhoq.*F2B)*ones(1,nknel)).*psi)].';
Bcontrib=[wf' *(((rhoq.*F2A.*jacobi)*ones(1,nknel)).*psi)].';

if elknrzb(1,3)==przb(3)
    CK=(rhoq.*Cfunk1a).' *wf;
else
    CK=0;
end

