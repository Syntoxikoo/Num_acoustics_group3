function []=checknormals(xyzb,topology,origo,varargin);
% []=checknormals(xyzb,topology,origo);
% function to check the local numbering in topology
% in order to produce consistent element normals

if nargin<2
   error('checknormals: insufficient number of input arguments');
end

if nargin==2
   origo=[0 0 0];
end

if nargin>3
   error('checknormals: too many input arguments');
end

IP=[0 0 0]; % Select point in center

[N, ncoltopo] = size(topology);
% I need to add the body number at a later stage
nknel=ncoltopo; % Number of nodes per element

Cdata=zeros(N,3);
rc=[];
for jj=1:N
   % Find the position and the normal vector at the center of the element
   elknxyzb=xyzb(topology(jj,1),:);
   for ikn=2:nknel
      elknxyzb=[elknxyzb; xyzb(topology(jj,ikn),:)];
   end
   
   [psi, xq, yq, zq, nx, ny, nz]=elemshape(elknxyzb,IP);
   % Are the radiusvector and the normal vector in the same hemi-space?
   rc=[rc ; xq yq zq]; % Save the center of the element for later use
   direction=sign((xq-origo(1))*nx+(yq-origo(2))*ny+(zq-origo(3))*nz);
   if direction < 0
      Cdata(jj,1)=1; % turn red if opposite direction
   else
      Cdata(jj,2)=1; % turn green if same direction
   end
end
figure;
xlabel('x');ylabel('y');zlabel('z');
patch('faces',topology,'vertices',xyzb,'FaceVertexCData',Cdata,'FaceColor','flat');
axis equal;
%for jj=1:N
%   text(rc(jj,1)*1.1,rc(jj,2)*1.1,rc(jj,3)*1.1,sprintf('%d',jj));
%end




