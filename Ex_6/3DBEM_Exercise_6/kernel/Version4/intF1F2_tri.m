function [F4A,F4B,CK]=intF1F2_tri(pxyzb,elknxyzb,k,singular)

%  [F2A,F2B,CK]=intF1F2_tri(pxyzb,elknxyzb,k,singular)
%  
%  Calculates the singular AND non-singular parts of the FA and FB
%  integrals.
%
%  Input:
%    pxyzb:    real vector containing the (x,y,z,body) values for
%              the point 'P'
%    elknxyzb: real matrix, one row for each node in the element
%              each row contains (x,y,z,body) for the node
%    k:        wavenumber real number
%    singular: If the node pxyzb is contained in the element singular 
%              denotes its position - otherwise singular is 0
%
%  Output:
%    F2A: complex vector, contains the contribution of each
%         node in the element to the B matrix (psi matrix)
%
%    F2B: complex vector, contains the contribution of each
%         node in the element to the A matrix (v matrix)
%  
%    CK:  real vector, contains the contribution of each
%         node in the element to the C constants

%  Based on my Ph.D. (pmj)

% Version for triangular elements - VCH 05-2007

% Singular integration extended to 6-node triangular elements - VCH 06-2012

persistent IP psi xq yq zq nx ny nz jacobi

nknel=size(elknxyzb,1);

%n=8; % Gaussrule order
n=16; % Gaussrule order
%n=40; % Gaussrule order
if isempty(IP)
   IPF=singgausstri(n,1:nknel);
end

[psi, xq, yq, zq, nx, ny, nz]=elemshapetri(elknxyzb,IPF(:,:,singular));
jacobi=sqrt(nx.^2+ny.^2+nz.^2);
%plot3(xq, yq, zq,'r.') % see the integration points on the 3D object plot

R=sqrt((xq-pxyzb(1)).^2+(yq-pxyzb(2)).^2+(zq-pxyzb(3)).^2);
Gfunk1=exp(-i*k*R)./R;
Bfunk1=(Gfunk1.*jacobi*ones(1,nknel)).*psi;
Cfunk1=-(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^2;
%Cfunk1=-(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^2.*(i*k*R+1).*Bfunk1;
Afunk1=((Cfunk1.*(i*k*R+1).*Gfunk1)*ones(1,nknel)).*psi;
F4A=IPF(:,3,singular)'*Afunk1;
F4B=IPF(:,3,singular)'*Bfunk1;

if elknxyzb(1,4)==pxyzb(4)
   CK=IPF(:,3,singular)'*(Cfunk1./R);
else
   CK=0;
end
