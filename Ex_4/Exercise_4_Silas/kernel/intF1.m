function [F1A,F1B,CK]=intF1(pxyb,elknxyb,k)

%  [F1A,F1B,CK]=intF1(pxyb,elknxyb,k)
%  
%  Calculates non-singular FA and FB integrals
%  for exterior collocation points 
%
%  Input variables:
%    pxyb:    real vector containing the (x,y,body) values for
%             the point 'P'.
%    elknxyb: real matrix, one row for each node in the element,
%             each row contains (x,y,body) for the node.
%    k: wavenumber.
%
%  Output:
%   -F1A: complex vector, contains the contribution of each
%         node in the element to the coefficient matrix for pressure.
%   -F1B: complex vector, contains the contribution of each 
%         node in the element to the coefficient matrix for normal velocity.
%   -CK:  contribution to the C constant.

% Vicente Cutanda Henriquez 5-2001.

nknel=size(elknxyb,1);

% Gauss-Legendre order
n=10;
[bp,wf]=gaussrule(n);

% obtain shape functions, normal vector and global coordinates of the integration points
[psi, xq, yq, nx, ny]=elemshape(elknxyb,bp);

jacobi=sqrt(nx.^2+ny.^2);    % jacobian = modulus of the normal vector

R1=sqrt((xq-pxyb(1)).^2+(yq-pxyb(2)).^2); % Distances from IPs to collocation point
dR1dn=((xq-pxyb(1)).*nx+(yq-pxyb(2)).*ny)./R1;


% Free field Green function and its derivative
G0dir=besselh(0,1,k*R1);
dG0dirdR1=besselh(1,1,k*R1);


%Calculates h1, h2 and h3
cons=-i*k*pi/2;
dG0dn=dG0dirdR1.*dR1dn;
h1=wf'*cons*(psi(:,1).*dG0dn);
h2=wf'*cons*(psi(:,2).*dG0dn);
h3=wf'*cons*(psi(:,3).*dG0dn);

% Finally we introduce the radiation term
Gtotal=(i/4)*G0dir;
g1P=wf'*(psi(:,1).*Gtotal.*jacobi);
g2P=wf'*(psi(:,2).*Gtotal.*jacobi);
g3P=wf'*(psi(:,3).*Gtotal.*jacobi);

F1A=[h1 h2 h3];
F1B=[g1P g2P g3P];


%we have to give dRdn1,R1,jacobi,wf
if pxyb(end)==elknxyb(1,end)
    CK=(dR1dn./R1)'*wf;
else
    CK=0;
end
