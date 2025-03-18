function p_field=FieldPnt3(fprzC,ps,vs,rzb,Topology,k,m,rho,c,varargin)

%  p_field=FieldPnt3(fprzC,ps,vs,rzb,Topology,k,m,rho,c,planeON);
%
%  Axisymetric BEM calculation.
%  Calculate pressure in field points given a solution on the boundary
%  
%  Input parameters (in SI units):
%  
%    -fprzC   : Field points' coordinates. Fisrt column is rho-coordinates,
%               second column is is z-coordinates and third column are the 
%               corresponding C constants (4*pi or 0)
%    -ps      : complex vector, the solution along the generator.
%    -vs      : normal velocity (complex column vector). Must be expanded
%               if the B matrix is expanded, see below.
%    -rzb     : Geometry matrix, one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 with the sign indicating
%               exterior/interior domain
%    -Topology: Topology matrix
%               One row for each element. Each row contains global node numbers
%               in the column no. corresponding to the local
%               node number. Last column contains the body number with the sign
%               indicating exterior/interior domain
%    -k       : Wave number (k=2*pi*f/c) (real scalar)
%    -m       : Circumferential mode number
%    -rho     : air density (Kg/m3). Zero value indicates that vs is a
%               velocity potential, instead of the normal velocity.
%    -c       : speed of sound (m/s)
%    -planeON : if 1, it assumes a reflecting plane at z=0, using Green's function
%
%  Output parameters:
%  
%    -p_field  : complex matrix, the sound pressure contribution of mode m on the field points.
%
% When non-continuous boundary conditions exist, the columns in
% matrix B can be expanded. An imaginary part must be added to
% the body numbers in 'topology' to indicate what boundary area
% the element belongs to (1,2,3...). The function 'bound2D' can be
% used to set up the equations and define the boundary areas.(VC)

% VCH 5-2016: A reflecting plane at z=0 can be defined.
% Could be improved by implementing the new Green's function into the
% integral functions, rather than calling them twice. An impedance term
% could also be included, like in 2DBEM.

% VCH 11-2018: B-matrix expansion, BC splitting. 


if nargin>9 && varargin{1}==1
   operate='y';
   if any(any(fprzC(:,2)<0))
       error('Field points cannot reach under z=0 when an impedance plane is defined');
   end
else
   operate='n';
end
[rzbnd,Topologynd,rzdum,topodum,rznodum]=nodummy(rzb,Topology,operate);
lndu=size(rzbnd,1);

% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
if isreal(Topology) | operate=='y' % calculations with impedance plane are not combined with B matrix extension
    BCtopo=Topologynd;
%     Bp=Ap;
else 
    BCelem=imag(Topology(:,end));
    Topologynd=real(Topology);
    Topology=real(Topology);
    % see help in function 'bound2D' (VC)
    [BCtopo,BCnodeA,BCnodeB]=bound2D(rzb,Topology,BCelem,'n');
end


M=size(rzb,1);
[Nfp,cc]=size(fprzC); 
if cc==2
   disp('Field points coordinates may include values of C constants (4*pi, domain, 0, off-domain)')
   disp('4*pi C constants are assumed')
   fprzC=[fprzC 4*pi*ones(Nfp,1)];
end
[nel, ncols]=size(Topology);
nknel=ncols-1;

p_field=zeros(Nfp,1);
for nfp=1:Nfp
   disp(['Field point row ' num2str(nfp) ' of ' num2str(Nfp)]);
   A=zeros(1,lndu);
   if exist('BCnodeB')
       B=zeros(1,size(BCnodeB,1));
   else
       B=zeros(1,lndu);
   end
   for iel=1:nel
      elknrzb(1:nknel,:)=rzb(Topology(iel,1:nknel),:);
      % Singular part
      [g,h,cjj] = intF2([fprzC(nfp,1:2) elknrzb(1,3)],elknrzb,m);
      if topodum(iel)~=0
          A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
          B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
      end
      if operate=='y'
          [g,h] = intF2([fprzC(nfp,1) -fprzC(nfp,2) elknrzb(1,3)],elknrzb,m);
          if topodum(iel)~=0
              A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
              B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
          end
      end
      % Oscillating part
      [g,h] = intF1([fprzC(nfp,1:2) elknrzb(1,3)],elknrzb,k,m);
      if topodum(iel)~=0
          A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
          B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
      end
      if operate=='y'
          [g,h] = intF1([fprzC(nfp,1) -fprzC(nfp,2) elknrzb(1,3)],elknrzb,k,m);
          if topodum(iel)~=0
              A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
              B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
          end
      end
   end
   if rho~=0, B= -i*k*rho*c*B; end % multiply B with -i*k*dens*c to get p instead of phi. Make rho zero to use phi.
   if fprzC(nfp,3)==0
      p_field(nfp,1) = NaN;
   else
      p_field(nfp,1) = (A*ps-B*vs)/fprzC(nfp,3); % compute field point
   end
end
