function [A,B,C]=point(nodesb,topologyb,k,xyz,varargin)

% [A,B,C]=point(nodesb,topologyb,k,xyz,nsingON,Tole)
%
% Calculates rows in A and B matrix form for the points xyz
% Usable for interior points (CHIEF) or exterior points
% If xyz contains several points, A and B are arrays, one row per point,
% and C is a vector, one value per point.
% This function checks if the points are too close to the boundary,
% and if they are near-singular integration is performed.
%
% Input variables:
%    -nodesb:    node x, y and z coordinates. One node per row.
%                Include body numbers.
%    -topologyb: numbers of the nodes in each element, ordered.
%                One element per row. The last column is the body number.
%                It is possible to include elements with different 
%                number of nodes (triangular or quadrilateral).
%                The remaining of the shorter rows should be set to NaN
%                and is nor read by the function.
%    -k:         wavenumber.
%    -xyz:       field or CHIEF points
%    -nsingON:   1, near-singular dealt with (default)
%                0, no near-singular check.
%    -Tole:      Allowed tolerance for near-singular check. It must be of 
%                the order of the expected minimum relative point-projection
%                distance. Default:1e-6. If Tole=0, no subdivision is made.
%
% Output variables:
%    -A:         rows of the coefficient matrix for the pressure.
%    -B:         rows of the coefficient matrix for the velocity.
%    -C:         C constants, solid angle as seen from the points.
%
% The arrays 'nodesb' and 'topology' must contain body numbers
% in the last column, with a sign to indicate
% exterior/interior (+/-) domain.(VC)
%
% When non-continuous boundary conditions exist, the columns in
% matrix B can be expanded. An imaginary part must be added to
% the body numbers in 'topology' to indicate what boundary area
% the element belongs to (0,1,2...). The function 'bound' can be
% used to set up the equations and define the boundary areas.(VC)

% Modified to include both triangular and quadrilateral elements (VCH 05-2007)

% Vicente Cutanda 12-2010: improvement of near-singular integrals and
% pre-check for increased speed and compatibility.

% Vicente Cutanda 2-2012: Input and help text modified.

tic;

if nargin==4
   nsingON=1;
   Tole=1e-6;
elseif nargin==5
   nsingON=varargin{1};
   Tole=1e-6;
elseif nargin==6
   nsingON=varargin{1};
   Tole=varargin{2};
elseif nargin>6
    error('Too many input arguments to function Point. Revise syntax.')
end


[M, ncolnodesb] = size(nodesb);
[np,cols]=size(xyz);
if ncolnodesb<4
   error('Error: The input arrays must include body numbers - Calculation aborted');
end
NumBodies = max(abs(nodesb(:,4)));

[N, nknel] = size(topologyb); % Limited to 3-node triangular and 4-node quadrilateral elements
nknel_quad=4;nknel_tri=3;

% "geometry": vector of length N (number of elements), where its values are 0 for a 
% quadrilateral element and 1 for a triangular element:
if nknel==5
    geometry=isnan(sum(topologyb,2));
    topotrib=topologyb;
elseif nknel==4
    geometry=ones(N,1);
    topotrib=[topologyb(:,1:end-1) NaN*ones(N,1) topologyb(:,end)];
else
    error('Element matrix input to Point is not well defined.')
end


%%%%%%%%%%%%%%%% A way to expand boundary conditions on triangular elements, deactivated (VCH)
Bv=0;

% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
% Not implemented for triangular elements
A=zeros(np,M);
if isreal(topologyb)
   BCtopo=topologyb;
   B=A;
else
   BCelem=imag(topologyb(:,end))+1;
   topologyb=real(topologyb);
   % see help in function 'bound' (VC)
   [BCtopo,BCnodeA,BCnodeB]=bound(nodesb,topologyb,BCelem,'n');
   B=zeros(np,size(BCnodeB,1));
end

if NumBodies ~= max(abs(topologyb(:,end)));
   error('Error: The input arrays must include consistent body numbers - Calculation aborted'); % ex/in (VC)
end

C=4*pi*ones(np,1);
for pp=1:np
   disp(sprintf('Field/CHIEF point calculation, point %g of %g',pp,np));
   for jj=1:N

      if geometry(jj)==0 %QUAD
         elknxyzb=nodesb(topologyb(jj,1:nknel_quad),:);
         xyz(pp,4)=elknxyzb(1,4); % this ensures body numbers always match and CKs are nonzero (VC)
         if nsingON==1
            is_close=nsingcheck(xyz(pp,:),elknxyzb);
         else
            is_close=0;
         end
         new=1; %New is true
         if ~is_close
            new=1; %New is true
            [F3A,F3B,CK]=intF3_quad(xyz(pp,:),elknxyzb,k,new);
            C(pp)=C(pp)+CK;
            A(pp,topologyb(jj,1:nknel_quad))=A(pp,topologyb(jj,1:nknel_quad))+F3A(1:nknel_quad);
            B(pp,BCtopo(jj,1:nknel_quad))=B(pp,BCtopo(jj,1:nknel_quad))+F3B(1:nknel_quad);
         else
            [F1A,F1B]=intF1_quad(xyz(pp,:),elknxyzb,k,new);
            [F4A,F4B,CK]=intF4_quad(xyz(pp,:),elknxyzb,Tole);
            C(pp)=C(pp)+CK;
            A(pp,topologyb(jj,1:nknel_quad))=A(pp,topologyb(jj,1:nknel_quad))+F1A(1:nknel_quad)+F4A(1:nknel_quad);
            B(pp,BCtopo(jj,1:nknel_quad))=B(pp,BCtopo(jj,1:nknel_quad))+F1B(1:nknel_quad)+F4B(1:nknel_quad);
         end
      end
      
      if geometry(jj)==1 %TRI
         elknxyzb=nodesb(topotrib(jj,1:nknel_tri),:);
         xyz(pp,4)=elknxyzb(1,4); % this ensures body numbers always match and CKs are nonzero (VC)
         if nsingON==1
            is_close=nsingcheck(xyz(pp,:),elknxyzb);
         else
            is_close=0;
         end
         if ~is_close
            new=1; %New is true
            [F3A,F3B,CK]=intF3_tri(xyz(pp,:),elknxyzb,k,new);
            C(pp)=C(pp)+CK;
            A(pp,topotrib(jj,1:nknel_tri))=A(pp,topotrib(jj,1:nknel_tri))+F3A(1:nknel_tri);
            if Bv
               B(pp)=B(pp)+F3B(1:nknel_tri)*v(jj,:)';
            else
               B(pp,topotrib(jj,1:nknel_tri))=B(pp,topotrib(jj,1:nknel_tri))+F3B(1:nknel_tri);
            end
         else
            [F4A,F4B,CK]=intF4_tri(xyz(pp,:),elknxyzb,k,Tole);
            C(pp)=C(pp)+CK;
            A(pp,topotrib(jj,1:nknel_tri))=A(pp,topotrib(jj,1:nknel_tri))+F4A(1:nknel_tri);
            if Bv
                B(pp)=B(pp)+F4B(1:nknel_tri)*v(jj,:)';
            else
                B(pp,BCtopo(jj,1:nknel_tri))=B(pp,BCtopo(jj,1:nknel_tri))+F4B(1:nknel_tri);
            end
         end

      end
   end
end

disp(['Calculation time for matrices of field points: ' num2str(toc)])

 
