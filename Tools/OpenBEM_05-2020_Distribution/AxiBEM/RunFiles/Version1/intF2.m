function [Bcontrib,Acontrib,CK]=intF2(przb,elknrzb,m)

%  [Bcontrib,Acontrib,CK]=intF2(przb,elknrzb,m)
%  
%  Calculates the singular part of the FA and FB integrals
%
%  Input:
%    przb:    real vector containing the (rho,z,body) values for
%             the point 'P'
%    elknrzb: real matrix, one row for each node in the element
%             each row contains (rho,z,body) for the node
%    m:       integer vector - circumferential mode number
%
%  Output:
%    Bcontrib: complex matrix, contains the contribution of each
%              node in the element to the B matrix. One column for
%              each m-value
%
%    Acontrib: complex matrix, contains the contribution of each
%              node in the element to the A matrix. One column for
%              each m-value.
%  
%    CK:       real vector, contains the contribution of each
%              node in the element to the C constants
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
%
% Inspiration for simultaneous m-values from:
%     Kuijpers, Verbeek & Verheij -
%     'An improved acoustic Fourier boundary element method
%      formulation using fast Fourier transform integration'

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

%  In my opinion the integration order needed for the generator integral
%  is determined by the order of shape function used when not dealing with (near)singular integrals (pmj)

[bp,wf]=rgauss(nknel-1,5);
x=bp;
[psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
% Note that 'nrho' and 'nz' have not been normalised, so
% The Jacobian should only be used for 'F1A' and 'F2A' (msj)
jacobi=sqrt(nrho.^2+nz.^2);

% Angular integration of singular part ('F2B' in thesis)
% Done analytically with elliptic integrals

% First see if integration needs to be refined due to near singularity
% or diagonal term
%threshold=0.2;
threshold=2;
dist2=(elknrzb(:,1)-przb(1)).^2+(elknrzb(:,2)-przb(2)).^2;
ell2=(elknrzb(1,1)-elknrzb(nknel,1)).^2+(elknrzb(1,2)-elknrzb(nknel,2)).^2;
[mindist2,closenode]=min(dist2);

if mindist2<ell2
   % Point 'P' is closer than one element length from the element
   % so it is necessary to check if the integration needs to be refined
   dist2=(rhoq-przb(1)).^2+(zq-przb(2)).^2; % Now dist2 contains the distances^2 from P to IPs
   [mindist2q,closeq]=min(dist2);

   if mindist2q<mindist2 % the closest point to P in the element is chosen, either a bp or a node
      mindist2=mindist2q; 
      singpos=x(closeq); % This should be improved as in 3DBEM; this "projection" is only approximate (VC)
   else
      singpos=closenode-2; % If the projection is on a node, singpos is exact (VC)
   end

   if mindist2<ell2*threshold^2 % a stricter check with the improved "mindist2"
      if mindist2<ell2*10*eps^2 % why not just mindist2<eps^2? (diagonal terms should give mindist2=0) (VC)
         % Diagonal term
         mindist2=1e-4*ell2; % a "fake" close distance is set to fool nsingrule; a truly singular integration should be implemented (VC)
      else
         % Near singularity
         disp('Near Singular Integration Performed');
      end

      [x,wf]=nsingrule(8,singpos,0.2*sqrt(mindist2/ell2));

      [psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x); % new IPs, new normals, etc.
      jacobi=sqrt(nrho.^2+nz.^2);

      if mindist2~=1e-4*ell2 & min((rhoq-przb(1)).^2+(zq-przb(2)).^2)<0.01*mindist2 % see 19 lines above:
                           % this check is needed to mend the error made
                           % there, since the new nsing IPs could be closer
                           % to P than the first rough ones, making the nsing ineffective (VC)
      % Yes it was a bit dirty of me to test for diagonal entries in this way (PMJ?) 
%         warning('Oblique collocation point too close to the boundary');
%         disp('PRESS RETURN TO CONTINUE');
%         pause;
      end
   end
end

% Angular integration
dz2=(zq-przb(2)).^2;
rstr2=(rhoq+przb(1)).^2+dz2;
kstr2=4./rstr2.*rhoq*przb(1);

forkval=1-(max(m)+1.01)^(-2); % Otherwise, Xm becomes negative due to round off
fork1=find(kstr2>forkval);
fork2=find(kstr2<=forkval);
F2A=zeros(length(zq),length(m));
F2B=zeros(length(zq),length(m));

if ~isempty(fork1)  % Use elliptic integrals
   % Deviations from the straight forward implementation of
   % Juhl's thesis are necessary due to finite precision arithmetic.
   % These are inspired by PMJ's Pascal file "h_gcos.pas" from the
   % FANPAC & DUCAT projects
   dz2=dz2(fork1);
   rstr2=rstr2(fork1);
   kstr2=kstr2(fork1);
   rstr=sqrt(rstr2);
   kstr=sqrt(kstr2);
   % kdash2=1-kstr2;
   % The following formula is more accurate though
   kdash2=((rhoq(fork1)-przb(1)).^2+dz2)./((rhoq(fork1)+przb(1)).^2+dz2);
   kdash=sqrt(kdash2);

   %  Elliptic integral (Matlab function, works on matrices):
   [kea,ea]=ellipke(kstr2);

   dkstrr=sqrt(przb(1)./rhoq(fork1))./rstr -kstr.*(rhoq(fork1)+przb(1))./rstr2;
   dkstrz=-kstr.*(zq(fork1)-przb(2))./rstr2;
   dke=(ea-kdash2.*kea)./kstr./kdash2;

   drstrr=(rhoq(fork1)+przb(1))./rstr;
   drstrz=(zq(fork1)-przb(2))./rstr;

   % Calculate the necessary In & dIn values
   mm=min(max(m),20); % series becomes unstable above this value
   for ii=mm:-1:0
      n=2*(mm-ii);
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
      if any(In<0) | any(In>pi/2 ./kdash)
         In, pi/2 ./kdash, n, mm
         error('In<0 or In>pi/(2*kdash)');
      end
      I(:,n+1)=In; % I(1)=I0 since matlab arrays must start at 1
      dI(:,n+1)=dIn;
   end

   if any(m>20) % Calculate Xm20,dXm20,Xm18,dXm18 for PMJ's extrapolation
      c=COEFFS(mm+1,1:2*mm+1);
      Xm20=I(:,2*mm+1);
      dXm20=dI(:,2*mm+1);
      for n=2*mm-2:-2:0
         Xm20=Xm20*4+c(n+1)*I(:,n+1);
         dXm20=dXm20*4+c(n+1)*dI(:,n+1);
      end
      Xm20=Xm20./2;
      dXm20=dXm20./2;
      mm=mm-2;
      c=COEFFS(mm+1,1:2*mm+1);
      Xm18=I(:,2*mm+1);
      dXm18=dI(:,2*mm+1);
      for n=2*mm-2:-2:0
         Xm18=Xm18*4+c(n+1)*I(:,n+1);
         dXm18=dXm18*4+c(n+1)*dI(:,n+1);
      end
      Xm18=Xm18./2;
      dXm18=dXm18./2;
      mm=mm+2;
   end

   % Calculate F2A and F2B
   for im=1:length(m)
      mm=m(im);
      if mm==0
         Xm=I(:,2*mm+1);
         dXm=dI(:,2*mm+1);
      elseif mm==18 & any(m>20)
         Xm=Xm18;
         dXm=dXm18;
      elseif mm==20 & any(m>20)
         Xm=Xm20;
         dXm=dXm20;
      elseif mm<=20
         c=COEFFS(mm+1,1:2*mm+1);
         Xm=I(:,2*mm+1);
         dXm=dI(:,2*mm+1);
         for n=2*mm-2:-2:0
            Xm=Xm*4+c(n+1)*I(:,n+1);
            dXm=dXm*4+c(n+1)*dI(:,n+1);
         end
         Xm=Xm./2;
         dXm=dXm./2;
      else % PMJ's extrapolation for high m-values
         signm=(-1)^mm;
         Xm=signm*Xm20.^((mm-18)/2).*Xm18.^((20-mm)/2).*(mm/20).^((mm-20)/180).*kstr2.^((20-mm));;
         dXm=signm*dXm20.^((mm-18)/2).*dXm18.^((20-mm)/2).*(mm/20).^((mm-20)/180).*kstr2.^((20-mm));;
%         Xm=signm*Xm20.^((mm-18)/2).*Xm18.^((20-mm)/2);
%         dXm=signm*dXm20.^((mm-18)/2).*dXm18.^((20-mm)/2);
      end

      % Test if Xm or dXm is illegal
      if any(abs(sign(Xm)-(-1)^mm)>1) | any(abs(sign(dXm)-(-1)^mm)>1)
         mm, rhoq(fork1), [Xm dXm]
         warning('Xm or dXm has the wrong sign!');
      end

      F2B(fork1,im)=(-1)^mm*4./(rstr2).*((rstr.*dXm.*dkstrr-Xm.*drstrr).*nrho(fork1)+(rstr.*dXm.*dkstrz-Xm.*drstrz).*nz(fork1));
      F2A(fork1,im)=((-1)^mm)*4./rstr.*Xm;
   end

   Cfunk1a(fork1,1)=4./(rstr2).*((rstr.*dke.*dkstrr-kea.*drstrr).*nrho(fork1)+(rstr.*dke.*dkstrz-kea.*drstrz).*nz(fork1));
end

if ~isempty(fork2) % Do the angular integration numerically

   nip=2*max(m)+8;
   bp_ang=[linspace(-1+1/nip,1-1/nip,nip)'];
   wf_ang=[2*ones(nip,1)./nip];
   if abs(sum(wf_ang)-2)>1e-4 | find(abs(bp_ang)>=1)
      sum(wf_ang), bp_ang
      error('Error in integration formula');
   end

   theta=bp_ang*(pi/2)+pi/2;
   bpones=ones(size(bp_ang));
   xones=ones(size(fork2));
   Ra=sqrt((rhoq(fork2).^2)*bpones'+przb(1)^2+((zq(fork2)-przb(2)).^2)*bpones' ...
      -2*przb(1)*rhoq(fork2)*cos(theta)');
   Bfr1=1./Ra;
   fr1=(-1./Ra.^3).* ...
       ((rhoq(fork2)*bpones'-przb(1)*xones*cos(theta)').*(nrho(fork2)*bpones')+ ...
       (zq(fork2)-przb(2)).*nz(fork2)*bpones');

   for im=1:length(m)
      mm=m(im);
      F2B(fork2,im)=2*fr1*(wf_ang.*cos(mm*theta))*pi/2;
      F2A(fork2,im)=2*Bfr1*(wf_ang.*cos(mm*theta))*pi/2;
   end
   Cfunk1a(fork2,1)=2*fr1*wf_ang*(pi/2);

%   F1B=fft([fr1 fr1(:,n_ang:-1:2)].').' *(2*pi/(2*n_ang));
%   F1A=fft([Bfr1 Bfr1(:,n_ang:-1:2)].').' *(2*pi/(2*n_ang));
%   Cfunk1=trapz(Cfr1.').' *2*pi/n_ang;
end


% Integration over element length
%Acontrib=[]; Bcontrib=[];
%for im=1:length(m)
%   Acontrib=[Acontrib [wf' *(((rhoq.*F2B(:,im))*ones(1,nknel)).*psi)].'];
%   Bcontrib=[Bcontrib [wf' *(((rhoq.*F2A(:,im).*jacobi)*ones(1,nknel)).*psi)].'];
%end
   
Acontrib=[]; Bcontrib=[];
for ikn=1:nknel
   Acontrib=[Acontrib; wf' *(((rhoq.*psi(:,ikn))*ones(1,length(m))).*F2B(:,1:length(m)))];
   Bcontrib=[Bcontrib; wf' *(((rhoq.*psi(:,ikn).*jacobi)*ones(1,length(m))).*F2A(:,1:length(m)))];
end

if elknrzb(1,3)==przb(3)
  CK=(rhoq.*Cfunk1a).' *wf;
else
  CK=0;
end

