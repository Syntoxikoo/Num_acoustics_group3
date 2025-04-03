function [A,B,CConst]=QTriEquat(xyzb,topologyb,k,varargin)

% [A,B,CConst]=QTriEquat(nodesb,topologyb,k,nsingON,Tole);
%
% Calculate C constants and A,B matrices for a given wavenumber k.
% This is the main function in the kernel. The matrix equation is:
%
%   0 = A*p - B*(-i*k*rho*c)*v + 4*pi*pI
%
% where p is the acoustic pressure, v is the particle velocity, rho is
% the density, c is the speed of sound and pI is the incident pressure.
% The exp(jwt) convention is used.
%
% Input variables:
%    -xyzb:      node x, y and z coordinates. One node per row.
%                Include body numbers.
%    -topologyb: numbers of the nodes in each element, ordered.
%                One element per row. The last column is the body number.
%                It is ONLY possible to include elements with 6
%                nodes (triangular).
%    -k:         wavenumber.
%    -nsingON:   1, near-singular integrals dealt with (default)
%                0, no near-singular check.
%    -Tole:      Minimum relative size of the subelements for near-singular
%                integrals. Default:1e-6. If Tole=0, no subdivision is made.
%
% Output variables:
%    -A:         coefficient matrix for the pressure. The C constants
%                are substracted from the diagonal.
%    -B:         coefficient matrix for the velocity.
%    -CConst:    C constants, solid angle as seen from the surface.
%
%
% The arrays 'xyzb' and 'topology' must contain body numbers
% in the last column, with a sign to indicate
% exterior/interior (+/-) domain.(VC)

% Peter M. Juhl 2000.

% Modified by Vicente Cutanda.

% Rene´ Christensen, 2002

% Vicente Cutanda 2-2006: generalization of near-singular integrals

% Vicente Cutanda 12-2010: improvement of near-singular integrals and
% pre-check for increased speed and compatibility.

% Vicente Cutanda 2-2012: Input and help text modified.

% Peter Juhl 5-2012: Quadratic triangular elements. Experimental
% formulation allowing only 6 node triangular elements

tic;

if nargin==3
    nsingON=1;
    Tole=1e-6;
elseif nargin==4
    nsingON=varargin{1};
    Tole=1e-6;
elseif nargin==5
    nsingON=varargin{1};
    Tole=varargin{2};
elseif nargin>5
    error('Too many input arguments to function TriQuadEquat. Revise syntax.')
end

[M, ncolxyzb] = size(xyzb);
if ncolxyzb<4
    error('Error: The input arrays must include body numbers - Calculation aborted'); % ex/in (VC)
end
NumBodies = max(abs(xyzb(:,4)));

[N, nknel] = size(topologyb); % Limited to 3-node ot 6-node triangular and 4-node quadrilateral elements
nknel_triq=6;

% "geometry": vector of length N (number of elements), where its values are 0 for a
% quadrilateral element and 1 for a triangular element:
if nknel==7
    topotrib=topologyb;
    Tolesing=1e-6;
else
    error('Element matrix input to QTriEquat is not well defined.')
end


%%%%%%%%%%%%%%%% A way to expand velocity boundary conditions on triangular elements, deactivated (VCH)
Bv=0;

% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
% Not implemented for triangular elements
A=zeros(M,M);
B=A;

if NumBodies ~= max(abs(topologyb(:,end)));
    error('Error: The input arrays must include consistent body numbers - Calculation aborted'); % ex/in (VC)
end

CConst=4*pi*(1+sign(xyzb(:,4)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
for jj=1:N
    disp(sprintf('BEM calculation, element %g of %g',jj,N));
    elknxyzb=xyzb(topotrib(jj,1:nknel_triq),:);
    
    % Deal with the singular nodes first
    singnodes=topotrib(jj,1:nknel_triq);
    % All singular integrals for the A matrix equals 0 for plane elements
    % This will not hold for the B matrix
    % We deal with singular nodes by subdivision - just like near singular
    % and using the same intF4_tri. Tolesing (see line 86) control the subdivision
    for ikn=1:nknel_triq
%        [F12A,F12B,CK]=intF4_tri(xyzb(singnodes(ikn),:),elknxyzb,k,Tolesing);
        [F12A,F12B,CK]=intF1F2_tri(xyzb(singnodes(ikn),:),elknxyzb,k,ikn);
        CConst(singnodes(ikn))=CConst(singnodes(ikn))+CK;
        A(singnodes(ikn),topotrib(jj,1:nknel_triq))=A(singnodes(ikn),topotrib(jj,1:nknel_triq))+F12A(1:nknel_triq);
        if Bv
            B(singnodes(ikn))=B(singnodes(ikn))+F12B(1:nknel_triq)*v(jj,:)';
        else
            B(singnodes(ikn),topotrib(jj,1:nknel_triq))=B(singnodes(ikn),topotrib(jj,1:nknel_triq))+F12B(1:nknel_triq);
        end
    end
    
    %Now deal with the non-singular nodes, and find near-singular nodes
    new=1; %New is true
    nsingdata=[];
    for ii=1:M
        %disp(sprintf('Node %g element %g',ii,jj));
        if isempty(find(singnodes==ii, 1)) % If ii is not among singnodes
            if nsingON==1
                is_close=nsingcheck(xyzb(ii,:),elknxyzb);
            else
                is_close=0;
            end
            if ~is_close
                % consider increasing integration points in intF3_tri
                [F3A,F3B,CK]=intF3_tri(xyzb(ii,:),elknxyzb,k,new);
                CConst(ii)=CConst(ii)+CK;
                A(ii,topotrib(jj,1:nknel_triq))=A(ii,topotrib(jj,1:nknel_triq))+F3A(1:nknel_triq);
                if Bv
                    B(ii)=B(ii)+F3B(1:nknel_triq)*v(jj,:)';
                else
                    B(ii,topotrib(jj,1:nknel_triq))=B(ii,topotrib(jj,1:nknel_triq))+F3B(1:nknel_triq);
                end
                new=0; %New is now false
            else
                nsingdata=[nsingdata ; ii];
            end
        end
    end
    
    % And finally deal with the near-singular nodes
    for ikn=1:size(nsingdata,1)
        ii=nsingdata(ikn);
        [F4A,F4B,CK]=intF4_tri(xyzb(ii,:),elknxyzb,k,Tole);
        CConst(ii)=CConst(ii)+CK;
        A(ii,topotrib(jj,1:nknel_triq))=A(ii,topotrib(jj,1:nknel_triq))+F4A(1:nknel_triq);
        if Bv
            B(ii)=B(ii)+F4B(1:nknel_triq)*v(jj,:)';
        else
            B(ii,topotrib(jj,1:nknel_triq))=B(ii,topotrib(jj,1:nknel_triq))+F4B(1:nknel_triq);
        end
    end
end

A=A-diag(CConst);

disp(['Calculation time for coefficient matrices: ' num2str(toc)])

