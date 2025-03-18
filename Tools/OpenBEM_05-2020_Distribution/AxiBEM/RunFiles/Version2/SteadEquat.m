function [A0,B0,CConst]=SteadEquat(rzb,Topology,beta,m)
%  Calculate C constants and A,B matrices:

[M, dummy]=size(rzb);
[nel, ncols]=size(Topology);
nknel=ncols-1;
A0=zeros(M,M);
B0=zeros(M,M);
h=zeros(nknel,1);
for ii=1:M
  cii=4*pi*(1+sign(rzb(ii,3)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
  disp(' ');
  disp(sprintf('A0 matrix calculation, row %g of %g',ii,M));
  for iel=1:nel
    elknrzb=rzb(Topology(iel,1),:);
    for ikn=2:nknel
       elknrzb=[elknrzb; rzb(Topology(iel,ikn),:)];
    end
    [g,h,cjj] = intF2(rzb(ii,:),elknrzb,m);
    cii=cii+cjj;
    for inode=1:nknel
       A0(ii,Topology(iel,inode))=A0(ii,Topology(iel,inode))+h(inode);
       B0(ii,Topology(iel,inode))=B0(ii,Topology(iel,inode))+g(inode);

    end
  end
  disp(sprintf('C constant %g = %1.15e',ii,cii/4/pi));
  CConst(ii)=cii/4/pi;
  A0(ii,ii)=A0(ii,ii)-cii;
end






