function p_field=FieldPnt3(fprzC,ps,vs,rzb,Topology,k,m,rho,c)

%  p_field=FieldPnt3(fprzC,ps,vs,rzb,Topology,k,m,rho,c);
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
%    -vs      : normal velocity (complex column vector)
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
%    -rho     : air density (Kg/m3)
%    -c       : speed of sound (m/s)
%
%  Output parameters:
%  
%    -p_field  : complex matrix, the sound pressure contribution of mode m on the field points.

M=size(rzb,1);
[Nfp,cc]=size(fprzC); 
if cc==2
   disp('Field points coordinates must include values of C constants (4*pi, domain, 0, off-domain)')
   disp('4*pi C constants are assumed')
   fprzC=[fprzC 4*pi*ones(Nfp,1)];
end
[nel, ncols]=size(Topology);
nknel=ncols-1;

p_field=zeros(Nfp,1);
for nfp=1:Nfp
   disp(['Field point row ' num2str(nfp) ' of ' num2str(Nfp)]);
   A=zeros(1,M);
   B=zeros(1,M);
   for iel=1:nel
      elknrzb(1:nknel,:)=rzb(Topology(iel,1:nknel),:);
      % Singular part
      [g,h,cjj] = intF2([fprzC(nfp,1:2) 0],elknrzb,m);
      A(Topology(iel,1:nknel))=A(Topology(iel,1:nknel))+h(1:nknel).';
      B(Topology(iel,1:nknel))=B(Topology(iel,1:nknel))+g(1:nknel).';
      % Oscillating part
      [g,h] = intF1([fprzC(nfp,1:2) 0],elknrzb,k,m);
      A(Topology(iel,1:nknel))=A(Topology(iel,1:nknel))+h(1:nknel).';
      B(Topology(iel,1:nknel))=B(Topology(iel,1:nknel))+g(1:nknel).';
   end
   B=i*k*rho*c*B; % multiply B with i*k*dens*c to get p instead of phi
   if fprzC(nfp,3)==0
      p_field(nfp,1) = NaN;
   else
      p_field(nfp,1) = (A*ps+B*vs)/fprzC(nfp,3); % compute field point
   end
end
