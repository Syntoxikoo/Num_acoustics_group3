function [F1A,F1B]=intF1_quad(pxyzb,elknxyzb,k,new)

%  [F1A,F1B]=intF1_quad(pxyzb,elknxyzb,k,new)
%  
%  Calculates the non-singular part of the FA and FB integrals
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
%    F1A: complex vector, contains the contribution of each
%         node in the element to the B matrix (psi matrix)
%
%    F1B: complex vector, contains the contribution of each
%         node in the element to the A matrix (v matrix)
%
%  Based on my Ph.D. (pmj)

persistent IP psi xq yq zq nx ny nz jacobi

nknel=size(elknxyzb,1);
%singular=0;         % assume that no singularity is present
n=8; % Gaussrule order
if isempty(IP)
   [bp,wf]=gaussrule(n);
   IP=[bp(ceil((1:n^2)/n)) bp(mod(0:n^2-1,n)+1) wf(ceil((1:n^2)/n)).*wf(mod(0:n^2-1,n)+1)];
end

if new
   [psi, xq, yq, zq, nx, ny, nz]=elemshape(elknxyzb,IP);
   jacobi=sqrt(nx.^2+ny.^2+nz.^2);
end
R=sqrt((xq-pxyzb(1)).^2+(yq-pxyzb(2)).^2+(zq-pxyzb(3)).^2);
Cfunk1=(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^3;
Cfunk1=Cfunk1.*(1-exp(-i*k*R).*(1+i*k*R));
Afunk1=(Cfunk1*ones(1,nknel)).*psi;
Bfunk1=(((exp(-i*k*R)-1).*jacobi./R)*ones(1,nknel)).*psi;
F1A=IP(:,3)'*Afunk1;
F1B=IP(:,3)'*Bfunk1;
