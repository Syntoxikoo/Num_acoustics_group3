function [nvect]=normals(xyzb,topologyb,see)

% [nvect]=normals(xyzb,topologyb,see)
%
% 3D version.
% Calculate normalvectors at the nodes using shape functions,
% element by element. Nodes belonging to two elements get the
% bisecting vector of their two normals.
% The vectors have unit length.
%
% Input: (rzb, topology) the geometry matrices with body numbers
%
% Output: (nvect) the normals, same rows as xyzb, 
%                 columns x, y and z coordinates
%
% If see='y', a figure is plotted with the normals as arrows (quiver).

% Vicente Cutanda Henriquez 04-2008


M = size(xyzb,1);
[N, ncols]=size(topologyb);
nknel=ncols-1;
if nknel==8
    x=[-1 -1;1 -1;1 1;-1 1;0 -1;1 0;0 1;-1 0];
elseif nknel==4
    x=[-1 -1;1 -1;1 1;-1 1];
elseif nknel==3
    x=[1 0;0 1;0 0];
elseif nknel==6
    x=[1 0;0 1;0 0;0.5 0.5;0 0.5;0.5 0];
end

nvect=zeros(M,3);

for iel=1:N
   if nknel==8 | nknel==4 
      [psi, xq, yq, zq, nx, ny, nz]=elemshape(xyzb(topologyb(iel,1:nknel),:),x);
   else
      [psi, xq, yq, zq, nx, ny, nz]=elemshapetri(xyzb(topologyb(iel,1:nknel),:),x);
   end
   nvect(topologyb(iel,1:nknel),:)=nvect(topologyb(iel,1:nknel),:) + [nx ny nz];
end

% normalize the vectors
nvect = nvect./(sqrt(nvect(:,1).^2 + nvect(:,2).^2 + nvect(:,3).^2)*[1 1 1]);

if see(1)=='y' | see(1)=='Y'
    plotsurf(xyzb,topologyb)
    hold on
    quiver3(xyzb(:,1),xyzb(:,2),xyzb(:,3),nvect(:,1),nvect(:,2),nvect(:,3));
end

