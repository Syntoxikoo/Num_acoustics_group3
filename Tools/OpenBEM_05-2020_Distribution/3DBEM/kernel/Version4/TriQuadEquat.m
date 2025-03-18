function [A,B,CConst]=TriQuadEquat(xyzb,topologyb,k,varargin)

% [A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON,Tole,PlaneON,Batch);
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
%    -planeON:   if 1, it assumes a reflecting plane at z=0, using Green's function
%    -Batch:     Disables (1) displayed progress for batch jobs (default 0).
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
% the element belongs to (1,2,3...). The function 'bound' can be
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

% VCH 3-2015: A and B matrices are filled using an intermediate variable
% Aele, Bele to increase speed in very large problems. Implemented for
% triangular elements only. Unused velocity expanxion code removed from the
% loop (TRI only).

% VCH 3-2016: B-matrix expansion brought back for triangular elements. The
% possibility of mixing quadrilateral and triangular elements is removed
% (it was driven by the control variable "geometry").

% CHANGES IN KERNEL FILES:
% nsing: line 32 in intF4tri.m  // check: line 22 in nsingcheck.m  // Ftest: line 125 in nsing2dTRIA.m  
% //  sing/nnsing: lines 37-48 in intF1F2tri.m  // gauss: line 30 in intF3_tri
% calcID='_nsing7_check2_Ftest4_sing16'; % >>> Large erors in the normal viscous velocity
%
% intF1f2_tri done with near-singular division instead of singular division, 
% and more IPs in non-singular integrals intF3_tri (VCH 7-2013)
% calcID='_nsing7_check2_Ftest4_nnsing7_gauss12'; 

% VCH 6-2016: A reflecting plane at z=0 can be defined (only triangular elements).
% Could be improved by implementing the new Green's function into the
% integral functions, rather than calling them twice. An impedance term
% could also be included, like in 2DBEM.


tic;

nsingON=1;
Tole=1e-6;
operate=0;
Batch=0; 
if nargin==4
   nsingON=varargin{1};
elseif nargin==5
   nsingON=varargin{1};
   Tole=varargin{2};
elseif nargin==6
   nsingON=varargin{1};
   Tole=varargin{2};
   operate=varargin{3};
elseif nargin==7
   nsingON=varargin{1};
   Tole=varargin{2};
   operate=varargin{3};
   Batch=varargin{4};
elseif nargin>7
    error('Too many input arguments to function TriQuadEquat. Revise syntax.')
end

if ~isreal(topologyb) & operate 
    error('B-matrix expansion is not implemented for the reflecting plane case. Only one of the two features can be used.');
end

[M, ncolxyzb] = size(xyzb);
if ncolxyzb<4
    error('Error: The input arrays must include body numbers - Calculation aborted'); % ex/in (VC)
end
NumBodies = max(abs(xyzb(:,4)));

[N, nknel] = size(topologyb); % Admits 3-node triangular, 6-node triangular and 4-node quadrilateral elements
nknel=nknel-1; % Should be either 4, 3 or 6 after this line;
if nknel~=4 && nknel~=3 && nknel~=6
    error('Error: The input arrays must include body numbers - Calculation aborted'); 
end

if NumBodies ~= max(abs(real(topologyb(:,end))));
    error('Error: The input arrays must include consistent body numbers - Calculation aborted'); % ex/in (VC)
end

if operate & nknel==4
    error('The reflecting plane case is only implemented for triangular elements.');
end


[xyzbnd,topologybnd,xyzdum,topodum,xyznodum]=nodummy3D(xyzb,topologyb,operate);
lndu=size(xyzbnd,1);

% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
A=zeros(lndu);
if isreal(topologyb) % && nknel~=5 
    BCtopo=topologybnd;
    B=A;
else % calculations with reflecting plane are not combined with B matrix extension
    BCelem=imag(topologyb(:,end));
    topologybnd=real(topologyb);
    topologyb=real(topologyb);
    % see help in function 'bound' (VC)
    [BCtopo,BCnodeA,BCnodeB]=bound(xyzb,topologyb,BCelem,'y');
    B=zeros(M,size(BCnodeB,1));
end

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CConst=4*pi*(1+sign(xyzbnd(:,4)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
for jj=1:N
     
    if ~Batch, disp(sprintf('BEM calculation, element %g of %g',jj,N)); end
    
    if nknel==4 %QUAD
        
        elknxyzb=xyzb(topologyb(jj,1:nknel),:);
        % Deal with the singular nodes first
        new=1; %New is true
        singnodes=topologyb(jj,1:nknel);
        
        for ikn=1:nknel
            [F1A,F1B]=intF1_quad(xyzb(singnodes(ikn),:),elknxyzb,k,new);
            [F2A,F2B,CK]=intF2_quad(xyzb(singnodes(ikn),:),elknxyzb,ikn);
            CConst(singnodes(ikn))=CConst(singnodes(ikn))+CK;    
            A(singnodes(ikn),topologyb(jj,1:nknel))=A(singnodes(ikn),topologyb(jj,1:nknel))+F1A(1:nknel)+F2A(1:nknel);
            B(singnodes(ikn),BCtopo(jj,1:nknel))=B(singnodes(ikn),BCtopo(jj,1:nknel))+F1B(1:nknel)+F2B(1:nknel);
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
                   A(ii,topologyb(jj,1:nknel))=A(ii,topologyb(jj,1:nknel))+F3A(1:nknel);
                   B(ii,BCtopo(jj,1:nknel))=B(ii,BCtopo(jj,1:nknel))+F3B(1:nknel);
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
            A(ii,topologyb(jj,1:nknel))=A(ii,topologyb(jj,1:nknel))+F1A(1:nknel)+F4A(1:nknel);
            B(ii,BCtopo(jj,1:nknel))=B(ii,BCtopo(jj,1:nknel))+F1B(1:nknel)+F4B(1:nknel);
            new=0; %New is now false
        end
    end
    
    
    if nknel==3 || nknel==6 % TRI3 and TRI6
        
        elenodes=topologyb(jj,1:nknel);
        elknxyzb=xyzb(elenodes,:);
        Aele=zeros(lndu,nknel);
        Bele=zeros(lndu,nknel); % Buffer matrices to increase speed
        
        % Deal with the singular nodes first
        new=1; %New is true
        % All singular integrals for the A matrix equals 0 for plane elements
        % This will not hold for the B matrix
        for ikn=1:nknel
            indx_nd=xyznodum(elenodes(ikn));
            if indx_nd~=0
                [F12A,F12B,CK12]=intF1F2_tri(xyzb(elenodes(ikn),:),elknxyzb,k,ikn);
                CConst(indx_nd)=CConst(indx_nd)+CK12*(1+xyzdum(indx_nd));
                if topodum(jj)~=0
                    Aele(indx_nd,:)=Aele(indx_nd,:)+F12A;
                    Bele(indx_nd,:)=Bele(indx_nd,:)+F12B;
                end
                if operate
                    if topodum(jj)~=0
                        [F12A,F12B]=intF1F2_tri([xyzb(elenodes(ikn),1:2) -xyzb(elenodes(ikn),3) xyzb(elenodes(ikn),4)],elknxyzb,k,ikn);
                        Aele(indx_nd,:)=Aele(indx_nd,:)+F12A;
                        Bele(indx_nd,:)=Bele(indx_nd,:)+F12B;
                    end
                end
            end
        end
        
        %Now deal with the non-singular nodes, and find near-singular nodes
        new=1; %New is true
        nsingdata=[];
        for ii=1:M
            %disp(sprintf('Node %g element %g',ii,jj));
            indx_nd=xyznodum(ii);
            if indx_nd~=0 && isempty(find(elenodes==ii, 1)) % If ii is not among elenodes
                if nsingON==1
                    is_close=nsingcheck(xyzb(ii,:),elknxyzb);
                else
                    is_close=0;
                end
                if ~is_close
                    [F3A,F3B,CK]=intF3_tri(xyzb(ii,:),elknxyzb,k,new); % consider increasing integration points in intF3_tri
                    CConst(indx_nd)=CConst(indx_nd)+CK*(1+xyzdum(indx_nd));
                    if topodum(jj)~=0
                        Aele(indx_nd,:)=Aele(indx_nd,:)+F3A;
                        Bele(indx_nd,:)=Bele(indx_nd,:)+F3B;
                        new=0; %New is now false
                    end
                    if operate
                        if topodum(jj)~=0
                            [F3A,F3B]=intF3_tri([xyzb(ii,1:2) -xyzb(ii,3) xyzb(ii,4)],elknxyzb,k,new); % consider increasing integration points in intF3_tri
                            Aele(indx_nd,:)=Aele(indx_nd,:)+F3A;
                            Bele(indx_nd,:)=Bele(indx_nd,:)+F3B;
                            new=0; %New is now false
                        end
                    end
                else
                    nsingdata=[nsingdata ; ii];
                end
            end
        end
        
        % And finally deal with the near-singular nodes
        for ikn=1:size(nsingdata,1)
            ii=nsingdata(ikn);
            indx_nd=xyznodum(ii);
            [F4A,F4B,CK]=intF4_tri(xyzb(ii,:),elknxyzb,k,Tole);
            CConst(indx_nd)=CConst(indx_nd)+CK*(1+xyzdum(indx_nd));
            if topodum(jj)~=0
                Aele(indx_nd,:)=Aele(indx_nd,:)+F4A;
                Bele(indx_nd,:)=Bele(indx_nd,:)+F4B;
            end
            if operate
                if topodum(jj)~=0
                    [F4A,F4B,]=intF4_tri([xyzb(ii,1:2) -xyzb(ii,3) xyzb(ii,4)],elknxyzb,k,Tole);
                    Aele(indx_nd,:)=Aele(indx_nd,:)+F4A;
                    Bele(indx_nd,:)=Bele(indx_nd,:)+F4B;
                end
            end
        end
        if topodum(jj)~=0
            A(:,topologybnd(topodum(jj),1:nknel))=A(:,topologybnd(topodum(jj),1:nknel))+Aele;
            %        B(:,elenodes)=B(:,elenodes)+Bele;
            B(:,BCtopo(topodum(jj),1:nknel))=B(:,BCtopo(topodum(jj),1:nknel))+Bele;
        end
        
    end
    
end

A=A-diag(CConst);
 
disp(['Calculation time for coefficient matrices: ' num2str(toc)])
 
