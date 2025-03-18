function [A,B]=HarmEquat(rzb,Topology,k,beta,m)
%  Calculate the frequency dependent part of the A and B matrices:

[M, dummy]=size(rzb);
[nel, ncols]=size(Topology);
nknel=ncols-1;
A=zeros(M,M);
B=zeros(M,M);
h=zeros(nknel,1);
for ii=1:M
  if M>=100
     disp(sprintf('A matrix calculation, row %3d of %3d',ii,M));
  end
  for iel=1:nel
    elknrzb=rzb(Topology(iel,1),:);
    for inode=2:nknel
       elknrzb=[elknrzb; rzb(Topology(iel,inode),:)];
    end
    [g,h] = intF1(rzb(ii,:),elknrzb,k,m);
    for inode=1:nknel
       A(ii,Topology(iel,inode))=A(ii,Topology(iel,inode))+h(inode);
       B(ii,Topology(iel,inode))=B(ii,Topology(iel,inode))+g(inode);
    end
  end
end






