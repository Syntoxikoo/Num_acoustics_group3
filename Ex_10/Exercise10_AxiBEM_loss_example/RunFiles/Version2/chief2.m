function [Achief,Bchief]=chief2(rzb,Topology,k,m,chiefpoint)

% [Achief,Bchief]=chief2(rzb,Topology,k,m,chiefpoint)
%
% Calculate an extra row for A,B matrices
% Helmholtz integral for a given CHIEF point in the interior
% The CHIEF method for remedying the non-uniqueness problem
%
% Input:
%   -rzb, Topology: Matrices difining geometry. Nodes and
%                   elements, with an extra column of body
%                   numbers.
%   -k:             Wavenumber.
%   -m:             Order of the circular expansion term.
%   -chiefpoint:    vector with rho, z and 1 (as body number),
%                   coordinates of the CHIEF point.


M=size(rzb,1);
[nel, ncols]=size(Topology);
nknel=ncols-1;
cii=4*pi;
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






