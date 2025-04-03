function [F1A,F1B]=intF1(przb,elknrzb,k,m)

%  [F1A,F1B]=intF1(przb,elknrzb,k,m)
%  
%  Calculates the non-singular part of the FA and FB integrals
%
%  Input:
%    przb:    real vector containing the (rho,z,body) values for
%             the point 'P'
%    elknrzb: real matrix, one row for each node in the element
%             each row contains (rho,z,body) for the node
%    k:       real scalar - the wave number
%    m:       integer scalar - circumferential mode number
%
%  Output:
%    F1A: complex vector, contains the contribution of each
%         node in the element to the B matrix
%
%    F1B: complex vector, contains the contribution of each
%         node in the element to the A matrix
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


[nknel,dummy]=size(elknrzb);


%%%%%%%%%%%%%%%%%%%%%%%%%% Nsing introduced here to improve calculation
%%%%%%%%%%%%%%%%%%%%%%%%%% with losses: Examine closer what it does! VCH 7.4.12
% Pre-check to see if the calculation point is close to the element, and
% obtain integration points accordingly
elknrz=elknrzb(:,1:2);prz(1,1:2)=przb(1:2);
maxnoddist=max(sqrt(sum(diff([elknrz; elknrz(1,:)]).^2,2))); % Get size of the maximum distance between nodes
maxdist=max(sqrt(sum((elknrz-ones(size(elknrz,1),1)*prz).^2,2))); % distance from the calculation point to the farthest node in the element:
if maxdist>maxnoddist*1.1 % it should not get nodes from adjoining elements *******************************
    [bp,wf]=gaussrule(30);
else
    
    % The function nsingrule takes care of near-singular AND singular (diagonal
    % coefficients) integrands, by means of a limited recursive interval
    % division. It may be convenient to implement a specialized singular
    % numerical integration rule to solve the singular case.
    [bp,wf]=nsingrule(8,elknrzb,przb);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     % Test plot
%     figure; plot(bp,zeros(length(bp)),'or');grid
%     title(['Node 1: ' num2str(elknrzb(1,1)) ' ; ' num2str(elknrzb(1,2)) ...
%            '  Node 3: ' num2str(elknrzb(3,1)) ' ; ' num2str(elknrzb(3,2)) ...
%            '  Point: ' num2str(przb(1)) ' ; ' num2str(przb(2))]);
%     xlabel(['bp points: ' num2str(length(bp))])
%     pause; close

%  In my opinion the integration order needed for the generator integral
%  is determined by the order of shape function used when not dealing with (near)singular integrals (pmj)

%%%%%%%%%%%%% [bp,wf]=gaussrule(30); % 10 >>> 30 %%%%%%%%%%%%%%%%%%%%
x=bp;
[psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
% Note that 'nrho' and 'nz' have not been normalised, so
% The Jacobian should only be used for 'F1A' and 'F2A' (msj)
jacobi=sqrt(nrho.^2+nz.^2);


% Angular integration of the non-singular part ('F1A' & 'F1B' in the thesis)
% Done numerically using base points and weights defined elsewhere.
if abs(k)<eps
   F1B=zeros(nknel,1);
   F1A=zeros(nknel,1);
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
   ang_int_fac=2;
   n_ang=max([ceil(abs(k)*dr_max*ang_int_fac); 200]);  % 10 >>> 100 >>>> 200 %%%%%%%%%%%%%%%%%%%%
   bp_ang=(-1+1/n_ang: 2/n_ang: 1-1/n_ang)';
   wf_ang=ones(n_ang,1)*(2/n_ang);

   a=0;
   b=pi;
   theta=(a+b)/2+(b-a)/2*bp_ang;
   bpones=ones(size(bp_ang));
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
   funk1=2*fr1*(wf_ang.*cos(m*theta))*(b-a)/2;
   Bfunk1=2*Bfr1*(wf_ang.*cos(m*theta))*(b-a)/2;
   F1B=rhoq.*funk1;
   F1A=rhoq.*Bfunk1.*jacobi;
    
   F1B=[wf' *((F1B*ones(1,nknel)).*psi)].';
   F1A=[wf' *((F1A*ones(1,nknel)).*psi)].';
end





