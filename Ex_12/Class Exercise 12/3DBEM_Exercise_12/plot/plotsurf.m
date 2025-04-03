function []=plotsurf(xyzb,topology,varargin);
%plotsurf(nodes,topology,elelight,nodered);
%Plots a surface defined by its nodes and topology matrix
%elelight (optional), has a list of indexes in topology. It then plots
%those elements in a lighter color.
%nodered (optional), has a list of indexes in xyzb. It then plots
%those nodes in red.

[N,nknel] = size(topology);
[M,cols]=size(xyzb);
if cols~=3
   xyz=xyzb(:,1:3);
   nknel=nknel-1;
else
   xyz=xyzb;
end

if nknel==8
   colstopo=[1 5 2 6 3 7 4 8];
else
   colstopo=1:nknel;
end
topologytmp=topology(:,colstopo);
nodestmp=xyzb(:,1:3);

Ncolor=ones(N,1)*[0.7 0.92 0.75];
if nargin>=3
   elelight=varargin{1};
   Ncolor(elelight)=Ncolor(elelight)*0.7;
end

figure;
patch('faces',topologytmp,'vertices',nodestmp,'FaceVertexCData',Ncolor,'FaceColor','flat');
hold on;
if nknel==8
   temp=topology(:,5:8);
   plot3(xyz(temp(1:end),1),xyz(temp(1:end),2),xyz(temp(1:end),3),'k.');
end
xlabel('x');ylabel('y');zlabel('z');

if nargin>=4
   nodered=varargin{2};
   plot3(xyz(nodered,1),xyz(nodered,2),xyz(nodered,3),'r.');
end

%axis equal;
rotate3d on;
grid