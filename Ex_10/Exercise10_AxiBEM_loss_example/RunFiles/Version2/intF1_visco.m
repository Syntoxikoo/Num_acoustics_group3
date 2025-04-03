function [F1A,F1B]=intF1_visco(przb,elknrzb,k,m,przn)

%  [F1A,F1B]=intF1(przb,elknrzb,k,m,przn)
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
%    przn:    normal vector at the collocation point (viscous mode calculation)
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

%  In my opinion the integration order needed for the generator integral
%  is determined by the order of shape function used when not dealing with (near)singular integrals (pmj)

[bp,wf]=gaussrule(10); % 10 >>> 30 %%%%%%%%%%%%%%%%%%%%
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
%     if (abs(imag(k)./real(k)))>0.1
%         % When the wavenumber has a significant imaginary part (viscous and
%         % thermal modes), the integral is calculated using the near-singular
%         % basepoints and weights centered at theta=0. (VCH 7-2011)
%         [bp_ang,wf_ang]=nsingrule(8,[-1 0 1;0 0 1;1 0 1],[-1 0 1]);
%         [bp,wf]=nsingrule(8,elknrzb,przb);x=bp;
%         [psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
%         jacobi=sqrt(nrho.^2+nz.^2);
%     else
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
        dr_max=max(sqrt((zq-przb(2)).^2+(rhoq+przb(1)).^2)-sqrt((zq-przb(2)).^2+(rhoq-przb(1)).^2));  %% ??? ((-) (+) (-) (-))
        ang_int_fac=2;
        n_ang=max([ceil(abs(k)*dr_max*ang_int_fac); 20]);  % 10 >>> 100 >>>> 200 %%%%%%%%%%%%%%%%%%%%
        bp_ang=(-1+1/n_ang: 2/n_ang: 1-1/n_ang)';
        wf_ang=ones(n_ang,1)*(2/n_ang);
%     end
    
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
    
    % Correction of vector directions of P and Q
    dotprodPQ=((nrho./jacobi)*przn(1))*cos(theta)' + ((nz./jacobi)*przn(2))*bpones'; 
%    fr1=fr1.*dotprodPQ;
    Bfr1=Bfr1.*dotprodPQ;
    
    funk1=2*fr1*(wf_ang.*cos(m*theta))*(b-a)/2;
    Bfunk1=2*Bfr1*(wf_ang.*cos(m*theta))*(b-a)/2;
    F1B=rhoq.*funk1;
    F1A=rhoq.*Bfunk1.*jacobi;
    
%     % Correction of vector directions of P and Q  >>>> theta dependence not
%     % yet implemented, to be done analytically?
     dotprodPQ=((nrho./jacobi)*przn(1)) + ((nz./jacobi)*przn(2));
%     F1A=F1A.*dotprodPQ;
     F1B=F1B.*dotprodPQ;
    
    F1B=[wf' *((F1B*ones(1,nknel)).*psi)].';
    F1A=[wf' *((F1A*ones(1,nknel)).*psi)].';
end

