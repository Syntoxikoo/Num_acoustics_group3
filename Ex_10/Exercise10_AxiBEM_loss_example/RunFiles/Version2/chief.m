function [A,B,chiefpoint]=chief(rzb,Topology,wavenum,beta,m)
% Calculate an extra row for A,B matrices
% Helmholtz integral for a random point in the interior
% The CHIEF method for remedying the non-uniqueness problem

k=wavenum;
[M, dummy]=size(rzb);
rmin=min(rzb(:,1));
rmax=max(rzb(:,1));
zmin=min(rzb(:,2));
zmax=max(rzb(:,2));
[nel, ncols]=size(Topology);
nknel=ncols-1;
cii=4*pi;
while cii>1e-3
  if rmin==0 & m==0
     chiefpoint=[0 zmin+rand*(zmax-zmin) 1];
  else
     chiefpoint=[rmin+rand*(rmax-rmin) zmin+rand*(zmax-zmin) 1];
  end
  A=zeros(1,M);
  B=zeros(1,M);
  h=zeros(nknel,1);
  cii=4*pi;
  disp(' ');
  disp(sprintf('CHIEF point calculation (%g,%g)',chiefpoint(1),chiefpoint(2)));
  for iel=1:nel
    elknrzb=rzb(Topology(iel,1),:);
    for ikn=2:nknel
       elknrzb=[elknrzb; rzb(Topology(iel,ikn),:)];
    end
    [g,h,cjj] = intF2(chiefpoint,elknrzb,m);
    cii=cii+cjj;
    for inode=1:nknel
       A(1,Topology(iel,inode))=A(1,Topology(iel,inode))+h(inode);
       B(1,Topology(iel,inode))=B(1,Topology(iel,inode))+g(inode);
    end
  end
  disp(sprintf('C constant = %1.15e',cii/4/pi));
end

% Calculate the non-singular part of the integral
for iel=1:nel
   elknrzb=rzb(Topology(iel,1),:);
   for ikn=2:nknel
      elknrzb=[elknrzb; rzb(Topology(iel,ikn),:)];
   end
   [g,h] = intF1(chiefpoint,elknrzb,k,m);
   for inode=1:nknel
      A(1,Topology(iel,inode))=A(1,Topology(iel,inode))+h(inode);
      B(1,Topology(iel,inode))=B(1,Topology(iel,inode))+g(inode);
   end
end






