function [Bcontrib,Acontrib]=intF1(przb,elknrzb,k,m)

%  [Bcontrib,Acontrib]=intF1(przb,elknrzb,k,m)
%  
%  Calculates the non-singular part of the FA and FB integrals
%
%  Input:
%    przb:    real vector containing the (rho,z,body) values for
%             the point 'P'
%    elknrzb: real matrix, one row for each node in the element
%             each row contains (rho,z,body) for the node
%    k:       real scalar - the wave number
%    m:       integer vector - circumferential mode numbers
%
%  Output:
%    Bcontrib: complex matrix. Each column contains the contribution
%              of the nodes in the element to the B matrix for an m-value
%
%    Acontrib: complex vector. Each column contains the contribution
%              of the nodes in the element to the A matrix for an m-value
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

nknel=size(elknrzb,1);

%  In my opinion the integration order needed for the generator integral
%  is determined by the order of shape function used when not dealing with (near)singular integrals (pmj)

[bp,wf]=rgauss(nknel-1,5);
%[bp,wf]=rgauss(nknel-1,20);
x=bp;
[psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
% Note that 'nrho' and 'nz' have not been normalised, so
% The Jacobian should only be used for 'F1A' and 'F2A' (msj)
jacobi=sqrt(nrho.^2+nz.^2);


% Angular integration of the non-singular part ('F1A' & 'F1B' in the thesis)
% Done numerically using base points and weights defined elsewhere.
if abs(k)<eps
   Acontrib=zeros(nknel,length(m));
   Bcontrib=zeros(nknel,length(m));
else
   % Hmm what order is needed for angular integration (pmj)
   % If we assume that the order needed is linear in ka_max (it's not)
   % we use ANGular_INTegration_FACtor to scale the angular integration
   % Rhoq determines the circular arc length
   % Note btw that oscillating integrands does not respond well
   % to very high order formulas (might as well subdivide)
   %
   % msj 991020:
   % What is important is max(k*(Rmax-Rmin)) rather than ka_max
   % This saves a lot of time at high frequencies
   % The mid-point formula is used, but should ideally be replaced
   % by a repeated 3 or 4 point Gauss formula
   dr_max=max(sqrt((zq-przb(2)).^2+(rhoq+przb(1)).^2) ...
             -sqrt((zq-przb(2)).^2+(rhoq-przb(1)).^2));
   ang_int_fac=2; % Double sampling frequency to avoid aliasing
   n_ang=max(2^ceil(log2((abs(k)*(dr_max+eps)+max(m))*ang_int_fac)),8);
   theta=pi*linspace(0,1,n_ang+1).';

   bpones=ones(size(theta));
   xones=ones(size(x));
   Ra=sqrt((rhoq.^2)*bpones'+przb(1)^2+((zq-przb(2)).^2)*bpones'...
      -2*przb(1)*rhoq*cos(theta)');
   fr1=((1-exp(-j*k*Ra).*(1+j*k*Ra))./Ra.^3).*...
       ((rhoq*bpones'-przb(1)*xones*cos(theta)').*(nrho*bpones')+...
       (zq-przb(2)).*nz*bpones');
%  Shorter way of writing the same expression
%   fr1=((1-exp(-j*k*Ra).*(1+j*k*Ra))./Ra.^3).*...
%       ((rhoq.*nrho+(zq-przb(2)).*nz)*bpones'-przb(1)*nrho*cos(theta)');
   Bfr1=(exp(-j*k*Ra)-1)./Ra;

   F1B=fft([fr1 fr1(:,n_ang:-1:2)].').' *(2*pi/(2*n_ang));
   F1A=fft([Bfr1 Bfr1(:,n_ang:-1:2)].').' *(2*pi/(2*n_ang));

   Acontrib=[]; Bcontrib=[];
   for ikn=1:nknel
      Acontrib=[Acontrib; wf' *(((rhoq.*psi(:,ikn))*ones(1,length(m))).*F1B(:,m+1))];
      Bcontrib=[Bcontrib; wf' *(((rhoq.*psi(:,ikn).*jacobi)*ones(1,length(m))).*F1A(:,m+1))];
   end    
end





