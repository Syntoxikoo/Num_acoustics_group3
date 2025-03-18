function [phi_field]=helmbemm_field(fprz,phi,dphidn,rzb,ElemNodeNum,k,m,varargin)

%  [phi_field]=helmbem(fprz,phi,dphidn,rzb,topology,k,m{,incoming})
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
%    -phi     : Complex column vector, the solution along the generator.
%    -dphidn  : The normal derivative of phi (i.e. normal velocity
%               at boundary).
%               A complex matrix with one row for each element
%               and one column for each local node in the element
%               A cell array of such matrices can be given when
%               there is more than one m-value
%    -k       : Wavenumber (k=2*pi*f/c) (real scalar)
%    -m       : Circumferential mode numbers (integer vector)
%    -rzb     : Geometry matrix
%               one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 (not necessary if there is
%               just one body)
%    -topology: Element topology (matrix of integers)
%               Empty for automatic generation
%    -incoming: Complex column vector where each row is the incoming
%               field at the corresponding point in 'fprz'
%               A cell array of such vectors can be given if there is
%               more than one m-value.
%               Can be empty if there is no incoming field
%
%  Output parameters:
%  
%    -phi_field: Complex column vector (if there is just one m-value)
%                Each row containing the solution in a field point.
%                If there is more than one m-value this is an array
%                of cells each containing such a column vector.

%  msj 000322
%  

%  Key points about input
[M, ncolrzb] = size(rzb);
if ncolrzb<3
	rzb=[rzb ones(M,1)];
end
NumBodies = rzb(M,3);

if isempty(ElemNodeNum)
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
        if inode>M
           ElemNodeNum(iel,3)=1;
           break;
        end
     end
   end
end
nel=size(ElemNodeNum,1);
nknel=size(ElemNodeNum,2)-1;

ik=1; % There is only one frequency at the moment

if ~iscell(phi)
   phi={phi};
   phi{ik,1:length(m)}=phi{1};
end

if ~iscell(dphidn)
   dphidn={dphidn};
   dphidn{ik,1:length(m)}=dphidn{1};
end
   

phi_field=FieldPnt(fprz,phi,dphidn,rzb,ElemNodeNum,k,m);

if nargin<8
   pI=[];
else
   pI=varargin{1};
end
if ~isempty(pI)
   if ~iscell(pI)
      pI={pI};
      pI{ik,1:length(m)}=pI{1};
   end
   for im=1:length(m)
      phi_field{ik,im}=phi_field{ik,im}+pI{ik,im};
   end
end

if length(m)==1 & length(k)==1
   phi_field=phi_field{1,1};
end



function [phi_field]=FieldPnt(fprz,phi,dphidn,rzb,Topology,k,m)
%  Calculate Phi in field points:

M=size(rzb,1);
nfpnt=size(fprz,1);
[nel, ncols]=size(Topology);
nknel=ncols-1;
cii=4*pi;
ik=1;

for ii=1:nfpnt
  Am=zeros(length(m),M);
  Bm=zeros(length(m),nel,nknel);
  disp(sprintf('Field point no. %g of %g',ii,nfpnt));
  for iel=1:nel
    for inode=1:nknel
       elknrzb(inode,:)=rzb(Topology(iel,inode),:);
    end
    % Singular part
    [g,h,cjj] = intF2([fprz(ii,:) 1],elknrzb,m);
%    cii=cii+cjj;
    for inode=1:nknel
       Am(:,Topology(iel,inode))=Am(:,Topology(iel,inode))+h(inode,:).';
       Bm(:,iel,inode)=g(inode,:).';
    end
    % Oscillating part
    [g,h] = intF1([fprz(ii,:) 1],elknrzb,k,m);
    for inode=1:nknel
       Am(:,Topology(iel,inode))=Am(:,Topology(iel,inode))+h(inode,:).';
       Bm(:,iel,inode)=Bm(:,iel,inode)+reshape(g(inode,:).',length(m),1,1);
    end
  end
  for im=1:length(m)
     phi_field{ik,im}(ii,1) = (Am(im,:)*phi{im}-sum(sum(squeeze(Bm(im,:,:)).*dphidn{im})))/cii;
  end
end








