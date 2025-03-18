function [F4A,F4B,CK]=intF4_tri(pxyzb,elknxyzb,k,Tole)

%  [F4A,F4B,CK]=intF4_tri(pxyzb,elknxyzb,k,Tole)
%  
%  Calculates non-singular FA and FB integrals
%  for very close points.
%
%  Input:
%    pxyzb:     real vector containing the (x,y,z,body) values for
%               the point 'P'
%    elknxyzb:  real matrix, one row for each node in the element
%               each row contains (x,y,z,body) for the node
%    k:         wavenumber real number
%    Tole:      Minimum relative size of the subelements for near-singular
%               integrals. Default:1e-6. If Tole=0, no subdivision is made.
%
%  Output:
%    F4A: complex vector, contains the contribution of each
%         node in the element to the B matrix (psi matrix)
%
%    F4B: complex vector, contains the contribution of each
%         node in the element to the A matrix (v matrix)
%
%  Based on my Ph.D. (pmj)

% Version for triangular elements - VCH 05-2007

% Vicente Cutanda 12-2010, near-singular integration improvements.

nknel=size(elknxyzb,1);
n=8; % Gaussrule order
ipp=find(sum((elknxyzb(:,1:3)-pxyzb(ones(1,nknel),1:3)).^2,2) <eps*100);
if length(ipp)~=1, error('Error finding singular node');end
IPF=singgausstri(n,ipp);

[psi, xq, yq, zq, nx, ny, nz]=elemshapetri(elknxyzb,IPF);
%plot3(xq, yq, zq,'r.') % see the integration points on the 3D object plot
jacobi=sqrt(nx.^2+ny.^2+nz.^2);

R=sqrt((xq-pxyzb(1)).^2+(yq-pxyzb(2)).^2+(zq-pxyzb(3)).^2);
Gfunk1=exp(-i*k*R)./R;
Bfunk1=(Gfunk1.*jacobi*ones(1,nknel)).*psi;
Cfunk1=-(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^2;
%Cfunk1=-(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^2.*(i*k*R+1).*Bfunk1;
Afunk1=((Cfunk1.*(i*k*R+1).*Gfunk1)*ones(1,nknel)).*psi;
F4A=IPF(:,3)'*Afunk1;
F4B=IPF(:,3)'*Bfunk1;

if elknxyzb(1,4)==pxyzb(4)
   CK=IPF(:,3)'*(Cfunk1./R);
else
   CK=0;
end
