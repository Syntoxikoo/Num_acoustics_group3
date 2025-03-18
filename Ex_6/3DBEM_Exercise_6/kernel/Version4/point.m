function [A,B,C]=point(nodesb,topologyb,k,xyz,varargin)

% [A,B,C]=point(nodesb,topologyb,k,xyz,nsingON,Tole,PlaneON)
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
%    -planeON    if 1, it assumes a reflecting plane at z=0, using Green's function
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

% VCH 3-2016: B-matrix expansion revised for triangular elements. The
% possibility of mixing quadrilateral and triangular elements is removed
% (it was driven by the control variable "geometry").

% VCH 6-2016: A reflecting plane at z=0 can be defined (only triangular elements).
% Could be improved by implementing the new Green's function into the
% integral functions, rather than calling them twice. An impedance term
% could also be included, like in 2DBEM.

tic;

nsingON=1;
Tole=1e-6;
operate=0;
if nargin==5
   nsingON=varargin{1};
elseif nargin==6
   nsingON=varargin{1};
   Tole=varargin{2};
elseif nargin==7
   nsingON=varargin{1};
   Tole=varargin{2};
   operate=varargin{3};
elseif nargin>7
    error('Too many input arguments to function TriQuadEquat. Revise syntax.')
end


if ~isreal(topologyb) && operate 
    error('B-matrix expansion is not implemented for the reflecting plane case. Only one of the two features can be used.');
end


[M, ncolnodesb] = size(nodesb);
[np,cols]=size(xyz);
if ncolnodesb<4
   error('Error: The input arrays must include body numbers - Calculation aborted');
end
NumBodies = max(abs(nodesb(:,4)));

[N, nknel] = size(topologyb); % Limited to 3-node triangular and 4-node quadrilateral elements
nknel=nknel-1; % Should be either 4, 3 or 6 after this line;
if nknel~=4 && nknel~=3 && nknel~=6
    error('Error: The input arrays must include body numbers - Calculation aborted'); 
end

if operate && nknel==4
    error('The reflecting plane case is only implemented for triangular elements.');
end

[xyzbnd,topologybnd,xyzdum,topodum,xyznodum]=nodummy3D(nodesb,topologyb,operate);
lndu=size(xyzbnd,1);


% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
A=zeros(np,lndu);
if isreal(topologyb)
   BCtopo=topologybnd;
   B=A;
else
   BCelem=imag(topologyb(:,end));
   topologybnd=real(topologyb);
   topologyb=real(topologyb);
   % see help in function 'bound' (VC)
   [BCtopo,BCnodeA,BCnodeB]=bound(nodesb,topologyb,BCelem,'n');
   B=zeros(np,size(BCnodeB,1));
end




if NumBodies ~= max(abs(topologyb(:,end)));
   error('Error: The input arrays must include consistent body numbers - Calculation aborted'); % ex/in (VC)
end

C=ones(np,1)*4*pi*(1+sign(nodesb(1,4)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
%C=4*pi*ones(np,1);

for pp=1:np
   disp(sprintf('Field/CHIEF point calculation, point %g of %g',pp,np));
   for jj=1:N

      if nknel==4 %QUAD
         elknxyzb=nodesb(topologyb(jj,1:nknel),:);
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
            A(pp,topologyb(jj,1:nknel))=A(pp,topologyb(jj,1:nknel))+F3A(1:nknel);
            B(pp,BCtopo(jj,1:nknel))=B(pp,BCtopo(jj,1:nknel))+F3B(1:nknel);
         else
            [F1A,F1B]=intF1_quad(xyz(pp,:),elknxyzb,k,new);
            [F4A,F4B,CK]=intF4_quad(xyz(pp,:),elknxyzb,Tole);
            C(pp)=C(pp)+CK;
            A(pp,topologyb(jj,1:nknel))=A(pp,topologyb(jj,1:nknel))+F1A(1:nknel)+F4A(1:nknel);
            B(pp,BCtopo(jj,1:nknel))=B(pp,BCtopo(jj,1:nknel))+F1B(1:nknel)+F4B(1:nknel);
         end
      end
      
      if nknel==3 || nknel==6 % TRI3 and TRI6
         elknxyzb=nodesb(topologyb(jj,1:nknel),:);
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
            if topodum(jj)~=0
                A(pp,topologybnd(topodum(jj),1:nknel))=A(pp,topologybnd(topodum(jj),1:nknel))+F3A(1:nknel);
                B(pp,BCtopo(topodum(jj),1:nknel))=B(pp,BCtopo(topodum(jj),1:nknel))+F3B(1:nknel);
            end
            if operate
                if topodum(jj)~=0
                    [F3A,F3B]=intF3_tri([xyz(pp,1:2) -xyz(pp,3) xyz(pp,4)],elknxyzb,k,new);
                    A(pp,topologybnd(topodum(jj),1:nknel))=A(pp,topologybnd(topodum(jj),1:nknel))+F3A(1:nknel);
                    B(pp,BCtopo(topodum(jj),1:nknel))=B(pp,BCtopo(topodum(jj),1:nknel))+F3B(1:nknel);
                end
            end
         else
            [F4A,F4B,CK]=intF4_tri(xyz(pp,:),elknxyzb,k,Tole);
            C(pp)=C(pp)+CK;
            if topodum(jj)~=0
                A(pp,topologybnd(topodum(jj),1:nknel))=A(pp,topologybnd(topodum(jj),1:nknel))+F4A(1:nknel);
                B(pp,BCtopo(topodum(jj),1:nknel))=B(pp,BCtopo(topodum(jj),1:nknel))+F4B(1:nknel);
            end
            if operate
                if topodum(jj)~=0
                    [F4A,F4B]=intF4_tri([xyz(pp,1:2) -xyz(pp,3) xyz(pp,4)],elknxyzb,k,Tole);
                    A(pp,topologybnd(topodum(jj),1:nknel))=A(pp,topologybnd(topodum(jj),1:nknel))+F4A(1:nknel);
                    B(pp,BCtopo(topodum(jj),1:nknel))=B(pp,BCtopo(topodum(jj),1:nknel))+F4B(1:nknel);
                end
            end
         end

      end
   end
end

disp(['Calculation time for matrices of field points: ' num2str(toc)])

 
