function [phi_field]=helmbem_field(fprz,phi,dphidn,k,rzb,incoming)

%  [phi_field]=helmbem(fprz,phi,dphidn,k,rzb,incoming)
%
%  Axisymetric BEM calculation.
%  Calculation of phi in random field points
%  given a solution on the boundary
%  
%  Input parameters (in SI units):
%  
%    -fprz    : Field point coordinates
%               first column is rho-coordinate
%               second column is z-coordinate
%    -phi     : complex vector, the solution along the generator.
%    -dphidn  : normal derivative of phi (complex column vector)
%               normal velocity of boundary = -1 * dphidn
%    -k       : Wavenumber (k=2*pi*f/c)
%    -rzb     : Geometry matrix
%               one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 (not necessary if there is
%               just one body)
%    -incoming: Function to calculate the incoming field from rzb
%               Usually an inline function
%               Alternatively the name of a function
%               [] or '' if there is no incoming field
%
%  Output parameters:
%  
%    -phi_field: complex vector, the solution in the field points.

%  msj 990811
%  

%  Key points about input
[M, ncolrzb] = size(rzb);
if ncolrzb<3
	rzb=[rzb ones(M,1)];
end
NumBodies = rzb(M,3);

inode=0;
iel=0;
for ibody=1:NumBodies
  inode=inode+1;
  while inode~=M & rzb(inode+1,3)==ibody
     iel=iel+1;
     if rzb(inode,:)==rzb(inode+1,:)
        inode=inode+1;
     end
     ElemNodeNum(iel,1)=inode;
     ElemNodeNum(iel,2)=inode+1;
     ElemNodeNum(iel,3)=inode+2;
     ElemNodeNum(iel,4)=ibody;
     inode=inode+2;
  end
end
N=iel;

% If 'incoming' is not an inline function, turn it into one
if isempty(incoming)
   incoming='';
end
if ~isa(incoming,'inline')
   if isempty(incoming) | incoming==''
      incoming=inline('wavenum*zeros(size(rzb,1),1)');
   else
      incoming=inline([incoming '(rzb,k)']);
   end
end
% Compute the incoming field in the nodes
pI=incoming(fprz,k);

% Circumferential mode number
m=0;
% Impedance (dummy at the moment)
beta=zeros(size(pI));

phi_field=FieldPnt(fprz,phi,dphidn,rzb,ElemNodeNum,k,beta,m);
phi_field=phi_field+pI;




function [phi_field]=FieldPnt(fprz,phi,dphidn,rzb,Topology,k,beta,m)
%  Calculate Phi in field points:

[M, dummy]=size(rzb);
[nfpnt, dummy]=size(fprz);
[nel, ncols]=size(Topology);
nknel=ncols-1;
cii=4*pi;

g=zeros(nknel,1);
h=zeros(nknel,1);
phi_field=zeros(nfpnt,1);
for ii=1:nfpnt
  A=zeros(1,M);
  B=zeros(1,M);
  disp(' ');
  disp(sprintf('Field point no. %g of %g',ii,nfpnt));
  for iel=1:nel
    elknrzb=rzb(Topology(iel,1),:);
    for inode=2:nknel
       elknrzb=[elknrzb; rzb(Topology(iel,inode),:)];
    end
    [g,h] = intFB([fprz(ii,1:2) 0],elknrzb,k,m);
    for inode=1:nknel
       A(Topology(iel,inode))=A(Topology(iel,inode))+h(inode);
       B(Topology(iel,inode))=B(Topology(iel,inode))+g(inode);
    end
  end
  phi_field(ii) = (A*phi-B*dphidn)/cii;
end








