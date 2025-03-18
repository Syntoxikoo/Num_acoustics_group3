function [A,B,CConst]=TriQuadEquat(xyzb,topologyb,k,varargin)

% [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON,Tole);
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
%                It is possible to include elements with different 
%                number of nodes (triangular or quadrilateral).
%                The remaining of the shorter rows should be set to NaN
%                and is nor read by the function. Triangular 6-node
%                elementes should be used over the whole geometry.
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
%
% When non-continuous boundary conditions exist, the columns in
% matrix B can be expanded. An imaginary part must be added to
% the body numbers in 'topology' to indicate what boundary area
% the element belongs to (0,1,2...). The function 'bound' can be
% used to set up the equations and define the boundary areas.(VC)

% Peter M. Juhl 2000.

% Modified by Vicente Cutanda.

% Rene´ Christensen, 2002

% Vicente Cutanda 2-2006: generalization of near-singular integrals

% Vicente Cutanda 12-2010: improvement of near-singular integrals and
% pre-check for increased speed and compatibility.

% Vicente Cutanda 2-2012: Input and help text modified.

% Peter Juhl 5-2012: Quadratic triangular elements. Experimental
% formulation allowing only 6 node triangular elements

% VCH 6-2012: quadratic triangular elements included in the regular formulation.

% CHANGES IN KERNEL FILES:
% nsing: line 32 in intF4tri.m  // check: line 22 in nsingcheck.m  // Ftest: line 125 in nsing2dTRIA.m  
% //  sing/nnsing: lines 37-48 in intF1F2tri.m  // gauss: line 30 in intF3_tri
% calcID='_nsing7_check2_Ftest4_sing16'; % >>> Large erors in the normal viscous velocity
%
% intF1f2_tri done with near-singular division instead of singular division, 
% and more IPs in non-singular integrals intF3_tri (VCH 7-2013)
% calcID='_nsing7_check2_Ftest4_nnsing7_gauss12'; 


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

[N, nknel] = size(topologyb); % Limited to 3-node triangular and 4-node quadrilateral elements
nknel_quad=4;%nknel_tri=3;nknel_triq=6;

% "geometry": vector of length N (number of elements), where its values are 0 for a 
% quadrilateral element and 1 for a triangular element:
if nknel==5
    geometry=isnan(sum(topologyb,2));
    topotrib=topologyb;
elseif nknel==4
    geometry=ones(N,1);
    topotrib=[topologyb(:,1:end-1) NaN*ones(N,1) topologyb(:,end)];
elseif nknel==7
    geometry=2*ones(N,1);
    topotrib=topologyb;
%    Tolesing=1e-6;
else
    error('Element matrix input to TriQuadEquat is not well defined.')
end


%%%%%%%%%%%%%%%% A way to expand velocity boundary conditions on triangular elements, deactivated (VCH)
Bv=0;



% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
% Not implemented for triangular elements
A=zeros(M,M);
if isreal(topologyb) && nknel~=5 
%    BCtopo=topologyb;
    B=A;
else
    BCelem=imag(topologyb(:,end))+1;
    topology=real(topologyb);
    % see help in function 'bound' (VC)
    [BCtopo,BCnodeA,BCnodeB]=bound(xyzb,topologyb,BCelem,'n');
    B=zeros(M,size(BCnodeB,1));
end

if NumBodies ~= max(abs(topologyb(:,end)));
    error('Error: The input arrays must include consistent body numbers - Calculation aborted'); % ex/in (VC)
end

CConst=4*pi*(1+sign(xyzb(:,4)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
for jj=1:N
     
    disp(sprintf('BEM calculation, element %g of %g',jj,N));
    
    if geometry(jj)==0 %QUAD
        
        elknxyzb=xyzb(topologyb(jj,1:nknel_quad),:);
        % Deal with the singular nodes first
        new=1; %New is true
        singnodes=topologyb(jj,1:nknel_quad);
        
        for ikn=1:nknel_quad
            [F1A,F1B]=intF1_quad(xyzb(singnodes(ikn),:),elknxyzb,k,new);
            [F2A,F2B,CK]=intF2_quad(xyzb(singnodes(ikn),:),elknxyzb,ikn);
            CConst(singnodes(ikn))=CConst(singnodes(ikn))+CK;    
            A(singnodes(ikn),topologyb(jj,1:nknel_quad))=A(singnodes(ikn),topologyb(jj,1:nknel_quad))+F1A(1:nknel_quad)+F2A(1:nknel_quad);
            B(singnodes(ikn),BCtopo(jj,1:nknel_quad))=B(singnodes(ikn),BCtopo(jj,1:nknel_quad))+F1B(1:nknel_quad)+F2B(1:nknel_quad);
            %new=0; %New is now false
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
                   [F3A,F3B,CK]=intF3_quad(xyzb(ii,:),elknxyzb,k,new);
                   CConst(ii)=CConst(ii)+CK;    
                   A(ii,topologyb(jj,1:nknel_quad))=A(ii,topologyb(jj,1:nknel_quad))+F3A(1:nknel_quad);
                   B(ii,BCtopo(jj,1:nknel_quad))=B(ii,BCtopo(jj,1:nknel_quad))+F3B(1:nknel_quad);
                   new=0; %New is now false
                else
                    nsingdata=[nsingdata ; ii];
                end
            end
        end
        
        % And finally deal with the near-singular nodes
        new=1; %New is true
        for ikn=1:size(nsingdata,1)
            ii=nsingdata(ikn);
            [F1A,F1B]=intF1_quad(xyzb(ii,:),elknxyzb,k,new);
            [F4A,F4B,CK]=intF4_quad(xyzb(ii,:),elknxyzb,Tole);
            CConst(ii)=CConst(ii)+CK;
            A(ii,topologyb(jj,1:nknel_quad))=A(ii,topologyb(jj,1:nknel_quad))+F1A(1:nknel_quad)+F4A(1:nknel_quad);
            B(ii,BCtopo(jj,1:nknel_quad))=B(ii,BCtopo(jj,1:nknel_quad))+F1B(1:nknel_quad)+F4B(1:nknel_quad);
            new=0; %New is now false
        end
    end
    
    
    if geometry(jj)==1 | geometry(jj)==2 % TRI3 and TRI6
        
        if geometry(jj)==1, nknel_tri=3; else nknel_tri=6; end 
        elknxyzb=xyzb(topotrib(jj,1:nknel_tri),:);
        
        % Deal with the singular nodes first
        new=1; %New is true
        singnodes=topotrib(jj,1:nknel_tri);
        % All singular integrals for the A matrix equals 0 for plane elements
        % This will not hold for the B matrix
        for ikn=1:nknel_tri
            [F12A,F12B,CK]=intF1F2_tri(xyzb(singnodes(ikn),:),elknxyzb,k,ikn); 
            CConst(singnodes(ikn))=CConst(singnodes(ikn))+CK;
            A(singnodes(ikn),topotrib(jj,1:nknel_tri))=A(singnodes(ikn),topotrib(jj,1:nknel_tri))+F12A(1:nknel_tri);
            if Bv
                B(singnodes(ikn))=B(singnodes(ikn))+F12B(1:nknel_tri)*v(jj,:)';
            else
                B(singnodes(ikn),topotrib(jj,1:nknel_tri))=B(singnodes(ikn),topotrib(jj,1:nknel_tri))+F12B(1:nknel_tri);
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
                   [F3A,F3B,CK]=intF3_tri(xyzb(ii,:),elknxyzb,k,new); % consider increasing integration points in intF3_tri
                   CConst(ii)=CConst(ii)+CK;    
                   A(ii,topotrib(jj,1:nknel_tri))=A(ii,topotrib(jj,1:nknel_tri))+F3A(1:nknel_tri);
                   if Bv
                       B(ii)=B(ii)+F3B(1:nknel_tri)*v(jj,:)';
                   else
                       B(ii,topotrib(jj,1:nknel_tri))=B(ii,topotrib(jj,1:nknel_tri))+F3B(1:nknel_tri);
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
            A(ii,topotrib(jj,1:nknel_tri))=A(ii,topotrib(jj,1:nknel_tri))+F4A(1:nknel_tri);
            if Bv
                B(ii)=B(ii)+F4B(1:nknel_tri)*v(jj,:)';
            else
                B(ii,topotrib(jj,1:nknel_tri))=B(ii,topotrib(jj,1:nknel_tri))+F4B(1:nknel_tri);
            end
        end
        
    end
    
end

A=A-diag(CConst);
 
disp(['Calculation time for coefficient matrices: ' num2str(toc)])
 
