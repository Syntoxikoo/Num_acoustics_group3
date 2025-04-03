function [Yavg]=AVGmem(rzbM)

%  [Yavg]=AVGmem(rzbM,topologyM)
%
%  Integrate displacement of a circular membrane using shape functions.
%
%  Input parameters:
%  
%    -rzbM    : cylindrical coordinates rho&z of the membrane nodes to integrate.
%               One row for each node, rho-coordinate in column 1,
%               z-coordinate in column 2. Body number in column 3 (2 and 3 not used)
%
%  Output parameters:
%
%    -Yavg :    row vector to apply to the displacement of the membrane nodes wi.
%               to obtain the mean displacement wm, a single number. 
%               wm = [Yavg] * [wi] 

%  Vicente Cutanda Henriquez, 10-2005


M=size(rzbM,1);

% Different kinds of elements can be defined on the supplied nodes
% by changing topology matrix:
%topologyM=[(1:3:(M-3))' (2:3:(M-2))' (3:3:(M-1))' (4:3:M)']; % cubic
%topologyM=[2*(1:(M-1)/2)'-1 2*(1:(M-1)/2)' 2*(1:(M-1)/2)'+1 ones((M-1)/2,1)]; % quadratic
topologyM=[(1:(M-1))' (2:M)' ones(M-1,1)]; % linear

% however, the Yavg vector obtained is different !!

[nel, ncols]=size(topologyM);
nknel=ncols-1;
[bp,wf]=gaussrule(20);

Yavg=zeros(1,M);
for iel=1:nel
    elknrzb(1:nknel,:)=rzbM(topologyM(iel,1:nknel),:);
	nknel=size(elknrzb,1);
	[psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,bp);
	jacobi=sqrt(nrho.^2+nz.^2);
	Ycontrib=[];
	for ikn=1:nknel
       Ycontrib=[Ycontrib; 2*pi*wf'*(rhoq.*psi(:,ikn).*jacobi)];
	end
   Yavg(1,topologyM(iel,1:nknel))=Yavg(1,topologyM(iel,1:nknel))+Ycontrib(1:nknel)';
end

