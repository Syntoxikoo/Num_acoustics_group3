function [p_field,fprB,fpzB,C]=FieldPnt(fpr,fpz,ps,vs,rzb,Topology,k,m,rho,c,varargin)

%  [p_field,fprB,fpzB,C]=FieldPnt(fpr,fpz,ps,vs,rzb,Topology,k,m,rho,c,ReplaceFP)
%
%  Axisymetric BEM calculation.
%  Calculate pressure in field points given a solution on the boundary
%  
%  Input parameters (in SI units):
%  
%    -fpr,fpz : Field points' coordinates, as a mesh. Two matrices of the
%               same size: fpr is rho-coordinates, fpz is z-coordinates
%    -ps      : complex vector, the solution along the generator.
%    -vs      : normal velocity (complex column vector)
%    -rzb     : Geometry matrix, one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 (not necessary if there is just one body)
%    -Topology: Topology matrix
%               One row for each element. Each row contains global node numbers
%               in the column no. corresponding to the local
%               node number. Last column contains the body number.
%    -k       : Wave number (k=2*pi*f/c) (real scalar)
%    -m       : Circumferential mode number
%    -rho     : air density (Kg/m3). Zero value indicates that vs is a
%               velocity potential, instead of the normal velocity.
%    -c       : speed of sound (m/s)
%    -ReplaceFP: If set to 1, checks for close surface points to replace field points.
%
%  Output parameters:
%  
%    -p_field  : complex matrix, the sound pressure contribution of mode m on the field points.
%    -fprB,fpzB: Field points' coordinates, as a mesh. Two matrices of the
%                same size: fpr is rho-coordinates, fpz is z-coordinates.
%                Corrected positions of on nodes close to the surface.
%    -C        : C constants from the calculation. 

if nargin<11, ReplaceFP=1; else, ReplaceFP=varargin{1}; end

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
   disp(['Field point row ' num2str(ff) ' of ' num2str(nfpff)]);
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
         A=zeros(1,M);
         B=zeros(1,M);
         cii=4*pi*(1+sign(rzb(1,3)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
         for iel=1:nel
            elknrzb(1:nknel,:)=rzb(Topology(iel,1:nknel),:);
            % Singular part
            [g,h,cjj] = intF2([fpr(ff,cc) fpz(ff,cc) elknrzb(1,3)],elknrzb,m);
            cii=cii+cjj;
            A(Topology(iel,1:nknel))=A(Topology(iel,1:nknel))+h(1:nknel).';
            B(Topology(iel,1:nknel))=B(Topology(iel,1:nknel))+g(1:nknel).';
            % Oscillating part
            [g,h] = intF1([fpr(ff,cc) fpz(ff,cc) elknrzb(1,3)],elknrzb,k,m);
            A(Topology(iel,1:nknel))=A(Topology(iel,1:nknel))+h(1:nknel).';
            B(Topology(iel,1:nknel))=B(Topology(iel,1:nknel))+g(1:nknel).';
         end
         C(ff,cc)=cii;
         if rho~=0, B= -i*k*rho*c*B; end % multiply B with -i*k*dens*c to get p instead of phi. Make rho zero to use phi.
         if abs(cii)<1e-3
             p_field(ff,cc) = 0;   % points outside the domain must have zero pressure
         else
             p_field(ff,cc) = (A*ps - B*vs)/cii; % compute field point
         end
      end
   end
end
