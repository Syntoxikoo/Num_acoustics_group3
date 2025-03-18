function [F2A,F2B,CK]=intF2_quad(pxyzb,elknxyzb,singular)

%  [F2A,F2B,CK]=intF2_quad(pxyzb,elknxyzb,singular)
%  
%  Calculates the singular part of the FA and FB integrals
%
%  Input:
%    pxyzb:    real vector containing the (x,y,z,body) values for
%             the point 'P'
%    elknxyzb: real matrix, one row for each node in the element
%             each row contains (x,y,z,body) for the node
%    singular: If the node pxyzb is contained in the element singular 
%             denotes its position - otherwise singular is 0
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
%
%  Based on my Ph.D. (pmj)

persistent IPs

nknel=size(elknxyzb,1);

% 'singular' is never 0, why then consider it?
n=8; % Gaussrule order
if isempty(IPs) | size(IPs,3)~=nknel
   IPs=singgauss2d(n,1:nknel);
end

[psi, xq, yq, zq, nx, ny, nz]=elemshape(elknxyzb,IPs(:,:,singular));
jacobi=sqrt(nx.^2+ny.^2+nz.^2);
   
R=sqrt((xq-pxyzb(1)).^2+(yq-pxyzb(2)).^2+(zq-pxyzb(3)).^2);
Cfunk1=-(nx.*(xq-pxyzb(1))+ny.*(yq-pxyzb(2))+nz.*(zq-pxyzb(3)))./R.^3;
Afunk1=(Cfunk1*ones(1,nknel)).*psi;
Bfunk1=((jacobi./R)*ones(1,nknel)).*psi;
CK=IPs(:,3,singular)'*Cfunk1;
F2A=IPs(:,3,singular)'*Afunk1;
F2B=IPs(:,3,singular)'*Bfunk1;
