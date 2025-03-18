function [F3A,F3B,CK]=intF3_quad(pxyzb,elknxyzb,k,new)

%  [F3A,F3B,CK,new]=intF3_quad(pxyzb,elknxyzb,k,new)
%  
%  Calculates non-singular FA and FB integrals
%  for exterior collocation points 
%
%  Input:
%    pxyzb:    real vector containing the (x,y,z,body) values for
%             the point 'P'
%    elknxyzb: real matrix, one row for each node in the element
%             each row contains (x,y,z,body) for the node
%	  k: wavenumber real number
%    new: if new=1 (true) then the element data is recalculated
%
%  Output:
%    F3A: complex vector, contains the contribution of each
%         node in the element to the B matrix (psi matrix)
%
%    F3B: complex vector, contains the contribution of each
%         node in the element to the A matrix (v matrix)
%
%  Based on my Ph.D. (pmj)

persistent IPF psi xq yq zq nx ny nz jacobi

nknel=size(elknxyzb,1);
%singular=0;         % assume that no singularity is present
n=8; % Gaussrule order

if isempty(IPF)
   [bp,wf]=gaussrule(n);
   IPF=[bp(ceil((1:n^2)/n)) bp(mod(0:n^2-1,n)+1) wf(ceil((1:n^2)/n)).*wf(mod(0:n^2-1,n)+1)];
end

if new
   [psi, xq, yq, zq, nx, ny, nz]=elemshape(elknxyzb,IPF);
   jacobi=sqrt(nx.^2+ny.^2+nz.^2);
end

R=sqrt((xq-pxyzb(1)).^2+(yq-pxyzb(2)).^2+(zq-pxyzb(3)).^2);
Gfunk1=exp(-i*k*R)./R;
Bfunk1=(Gfunk1.*jacobi*ones(1,nknel)).*psi;
Cfunk1=-(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^2;
%Cfunk1=-(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^2.*(i*k*R+1).*Bfunk1;
Afunk1=((Cfunk1.*(i*k*R+1).*Gfunk1)*ones(1,nknel)).*psi;
F3A=IPF(:,3)'*Afunk1;
F3B=IPF(:,3)'*Bfunk1;

if elknxyzb(1,4)==pxyzb(4)
   CK=IPF(:,3)'*(Cfunk1./R);
else
   CK=0;
end
