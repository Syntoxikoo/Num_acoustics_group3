function [p_field,fprB,fpzB,C]=FieldPnt(fpr,fpz,ps,vs,rzb,Topology,k,m,rho,c,varargin)

%  [p_field,fprB,fpzB,C]=FieldPnt(fpr,fpz,ps,vs,rzb,Topology,k,m,rho,c,ReplaceFP,planeON)
%
%  Axisymetric BEM calculation.
%  Calculate pressure in field points given a solution on the boundary
%  
%  Input parameters (in SI units):
%  
%    -fpr,fpz : Field points' coordinates, as a mesh. Two matrices of the
%               same size: fpr is rho-coordinates, fpz is z-coordinates
%    -ps      : complex vector, the solution along the generator.
%    -vs      : normal velocity (complex column vector). Must be expanded
%               if the B matrix is expanded, see below.
%    -rzb     : Geometry matrix, one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 (not necessary if there is just one body)
%    -Topology: Topology matrix
%               One row for each element. Each row contains global node numbers
%               in the column no. corresponding to the local
%               node number. Last column contains the body number.
%    -k       : Wave number (k=2*pi*f/c) (real scalar). 
%    -m       : Circumferential mode number
%    -rho     : air density (Kg/m3). Zero value indicates that vs is a
%               velocity potential, instead of the normal velocity.
%    -c       : speed of sound (m/s)
%    -ReplaceFP: If set to 1, checks for close surface points to replace field points.
%    -planeON : if 1, it assumes a reflecting plane at z=0, using Green's function
%
%  Output parameters:
%  
%    -p_field  : complex matrix, the sound pressure contribution of mode m on the field points.
%    -fprB,fpzB: Field points' coordinates, as a mesh. Two matrices of the
%                same size: fpr is rho-coordinates, fpz is z-coordinates.
%                Corrected positions of on nodes close to the surface.
%    -C        : C constants from the calculation. 
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

if nargin<11, ReplaceFP=1; else, ReplaceFP=varargin{1}; end

if nargin>11 && varargin{2}==1
   operate='y';
   if any(any(fpz<0))
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
[nfpff,nfpcc]=size(fpr); 
if ~all(size(fpr)==size(fpz))
   error('Field point matrices of coordinates mismatch')
end
[nel, ncols]=size(Topology);
nknel=ncols-1;
fprB=fpr;fpzB=fpz;
cii=0;
g=zeros(nknel,1);
h=zeros(nknel,1);
p_field=zeros(size(fpz));
C=zeros(size(fpz));
for ff=1:nfpff
   for cc=1:nfpcc
      fprT=fpr;fprT(ff,cc)=NaN;fpzT=fpz;fpzT(ff,cc)=NaN;
      minFPdist=min(min(sqrt((fpr(ff,cc)-fprT).^2 + (fpz(ff,cc)-fpzT).^2)));
      Spoint=find(sqrt((fpr(ff,cc)-rzb(:,1)).^2 + (fpz(ff,cc)-rzb(:,2)).^2)<0.3*minFPdist);
      if ~isempty(Spoint)& ReplaceFP
         [dd,Spoint1]=min(sqrt((fpr(ff,cc)-rzb(:,1)).^2 + (fpz(ff,cc)-rzb(:,2)).^2));
         p_field(ff,cc) = ps(Spoint1);
         fprB(ff,cc)=rzb(Spoint1,1);fpzB(ff,cc)=rzb(Spoint1,2);
         disp('Field point replaced by surface point')
         disp(['Distance, field point to nearest field point: ' num2str(minFPdist*1000) ' mm'])
         disp(['Distance, field point to replacing surface point: ' ...
            num2str(sqrt((fpr(ff,cc)-rzb(Spoint1,1)).^2 + (fpz(ff,cc)-rzb(Spoint1,2)).^2)*1000) ' mm'])
      else
          A=zeros(1,lndu);
          if exist('BCnodeB')
              B=zeros(1,size(BCnodeB,1));
          else
              B=zeros(1,lndu);
          end
         cii=4*pi*(1+sign(rzb(1,3)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
         for iel=1:nel
            elknrzb(1:nknel,:)=rzb(Topology(iel,1:nknel),:);
            % Singular part
            [g,h,cjj] = intF2([fpr(ff,cc) fpz(ff,cc) elknrzb(1,3)],elknrzb,m);
            cii=cii+cjj;
            if topodum(iel)~=0
                A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
                B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
            end
            if operate=='y'
                [g,h] = intF2([fpr(ff,cc) -fpz(ff,cc) elknrzb(1,3)],elknrzb,m);
                if topodum(iel)~=0
                    A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
                    B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
                end
            end
            % Oscillating part
            [g,h] = intF1([fpr(ff,cc) fpz(ff,cc) elknrzb(1,3)],elknrzb,k,m);
            if topodum(iel)~=0
                A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
                B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
            end
            if operate=='y'
                [g,h] = intF1([fpr(ff,cc) -fpz(ff,cc) elknrzb(1,3)],elknrzb,k,m);
                if topodum(iel)~=0
                    A(Topologynd(topodum(iel),1:nknel))=A(Topologynd(topodum(iel),1:nknel))+h(1:nknel).';
                    B(BCtopo(topodum(iel),1:nknel))=B(BCtopo(topodum(iel),1:nknel))+g(1:nknel).';
                end
            end
         end
         C(ff,cc)=cii;
         if rho~=0, B= -1i*k*rho*c*B; end % multiply B with -i*k*dens*c to get p instead of phi. Make rho zero to use phi.
         if abs(cii)<1e-3
             p_field(ff,cc) = 0;   % points outside the domain must have zero pressure
         else
             p_field(ff,cc) = (A*ps - B*vs)/cii; % compute field point
         end
      end
   end
   disp(['Field point row ' num2str(ff) ' of ' num2str(nfpff) '   Last C constant: ' num2str(cii)]);
end
