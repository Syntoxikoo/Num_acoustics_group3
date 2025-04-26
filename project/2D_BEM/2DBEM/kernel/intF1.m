function [F1A,F1B,CK]=intF1(pxyb,elknxyb,k,betaP)

%  [F1A,F1B,CK]=intF1(pxyb,elknxyb,k,betaP)
%  
%  Calculates non-singular FA and FB integrals
%  for exterior collocation points 
%
%  Input variables:
%    pxyb:    real vector containing the (x,y,body) values for
%             the point 'P'.
%    elknxyb: real matrix, one row for each node in the element,
%             each row contains (x,y,body) for the node.
%	  k: wavenumber.
%   -betaP:   admittance of the impedance plane. If its value is 0,
%             a rigid plane is considered (infinite impedance); if its value
%             is NaN (not-a-number), free-field is assumed.
%
%  Output:
%   -F1A: complex vector, contains the contribution of each
%         node in the element to the coefficient matrix for pressure.
%   -F1B: complex vector, contains the contribution of each 
%         node in the element to the coefficient matrix for normal velocity.
%   -CK:  contribution to the C constant.

% Vicente Cutanda 5-2001.

% Peter Juhl, 3-2002, added singular integration of dia.terms using Vicente's 
% near singular integration scheme.

persistent bpregular wfregular bpleft wfleft bpmiddle wfmiddle bpright wfright

nknel=size(elknxyb,1);

% Gauss-Legendre order
n=10;

if isempty(bpregular)
   [bpregular,wfregular]=gaussrule(n);
   [bpleft,wfleft]=nsingrule(8,-1,0.002);
   [bpmiddle,wfmiddle]=nsingrule(8,0,0.002);
   [bpright,wfright]=nsingrule(8,1,0.002);
end

% check for diagonal entry
nknel=size(elknxyb,1);
pntmat=ones(nknel,1)*pxyb;
pnttest=sum(abs((elknxyb-pntmat)'));
whatpnt=find(pnttest<100*eps);
if isempty(whatpnt) 
    bp=bpregular;
    wf=wfregular;
elseif whatpnt==1
    bp=bpleft;
    wf=wfleft;
elseif whatpnt==nknel
    bp=bpright;
    wf=wfright;
else
    bp=bpmiddle;
    wf=wfmiddle;
end

% NEAR-SINGULAR INTEGRALS (VCH 10-2014)
% Check to see if the calculation point is close to the element, and
% obtain integration points accordingly
elknxy=elknxyb(:,1:2);pxy(1,1:2)=pxyb(1:2);
maxnoddist=max(sqrt(sum(diff([elknxy; elknxy(1,:)]).^2,2))); 
% Get size of the maximum distance between nodes
maxdist=max(sqrt(sum((elknxy-ones(size(elknxy,1),1)*pxy).^2,2))); 
% Distance from the calculation point to the farthest node in the element:
if maxdist<maxnoddist*1.1 % it should not get nodes from adjoining elements
    % The function nsingrule takes care of near-singular integrands, by 
    % means of a limited recursive interval division.
    [bp,wf]=nsingruleP(8,elknxyb,pxyb);
end



% obtain shape functions, normal vector and global coordinates of the integration points
[psi, xq, yq, nx, ny]=elemshape(elknxyb,bp);

jacobi=sqrt(nx.^2+ny.^2);    % jacobian = modulus of the normal vector
%nx=nx./jacobi; ny=ny./jacobi; % normalize normal vector

R1=sqrt((xq-pxyb(1)).^2+(yq-pxyb(2)).^2); % Distances from IPs to collocation point
dR1dn=((xq-pxyb(1)).*nx+(yq-pxyb(2)).*ny)./R1;

if sum(abs(elknxyb(:,2)))~=0 | isnan(betaP)  % avoid dummy elements, which have all y = 0

   R2=sqrt((xq-pxyb(1)).^2+(yq+pxyb(2)).^2); % Distances from IPs to image collocation point
   dR2dn=((xq-pxyb(1)).*nx+(yq+pxyb(2)).*ny)./R2;

   [G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k,pxyb,betaP,xq,yq);

   %Calculates h1, h2 and h3 without the corrective term
   cons=-i*k*pi/2;
   dG0dn=dG0dirdR1.*dR1dn+dG0refdR2.*dR2dn;
   h1=wf'*cons*(psi(:,1).*dG0dn);
   h2=wf'*cons*(psi(:,2).*dG0dn);
   h3=wf'*cons*(psi(:,3).*dG0dn);

   %Now we have to use the derivatives of the correction term
   %to complete h1, h2, h3
   dPbetadn=dPbetadx.*nx + dPbetady.*ny;
   h1P=2*pi*wf'*(psi(:,1).*dPbetadn);
   h2P=2*pi*wf'*(psi(:,2).*dPbetadn);
   h3P=2*pi*wf'*(psi(:,3).*dPbetadn);

   % Finally we introduce the radiation term
   Gtotal=(i*pi/2)*(G0dir+G0ref)+(2*pi*Pbeta); % Bug, VCH 7-2011
%   Gtotal=(i/4)*(G0dir+G0ref)+Pbeta;
   g1P=wf'*(psi(:,1).*Gtotal.*jacobi);
   g2P=wf'*(psi(:,2).*Gtotal.*jacobi);
   g3P=wf'*(psi(:,3).*Gtotal.*jacobi);

   F1A=[h1+h1P h2+h2P h3+h3P];
   F1B=[g1P g2P g3P];

else
   F1A=zeros(1,size(elknxyb,1));
   F1B=zeros(1,size(elknxyb,1));
end

%we have to give dRdn1,R1,jacobi,wf
if pxyb(end)==elknxyb(1,end)
   CK=(dR1dn./R1)'*wf;
else
   CK=0;
end
