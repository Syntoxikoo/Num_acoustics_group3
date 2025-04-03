function []=playresult(xyzb,topology,result,varargin)
% []=playresult(xyzb,topology,result,fps)
% Plays a movie of result as color on a patch surface with changing phase.
%
% xyzb : nodal coordinates
% topology: relation between elements and nodes
% result : complex vector to be plotted
% fps : frames per second (default: 5)

% VC 11-2000

[N,nknel] = size(topology);
[M,cols]=size(xyzb);
if cols~=3
   xyz=xyzb(:,1:3);
   nknel=nknel-1;
else
   xyz=xyzb;
end

if nargin>3
   fps=varargin{1};
else
   fps=5;
end

disp(' Creating the frames...');
figure;
res=abs(result);
drawone(xyz,topology,res,nknel)
[camin camax]=caxis;
caxis([-camax camax]);
colorbar vert;
xlabel('x');ylabel('y');zlabel('z');
%axis equal;
%rotate3d on;
view(-25,40);

reso=20;
for j = 1:reso
   res=abs(result).*cos(angle(result)+2*pi/reso*j);
   drawone(xyz,topology,res,nknel)
   F(j) = getframe;
end
disp(' Click on the figure to start the movie');
[xx,yy,bb]=ginput(1);

% Play the movie
bb=1;
while bb==1
   movie(F,1,fps); 
   disp(' Left click: once more,  right click: end and close');
   [xx,yy,bb]=ginput(1);
end

close;


function []=drawone(xyz,topology,result,nknel)

if nknel==8
   [dummy,ii]=min([abs(abs(topology(:,7))-abs(topology(:,5))) abs(abs(topology(:,6))-abs(topology(:,8)))]');
   ind1=find(ii==1);
   ind2=find(ii==2);
   topology(ind1,[10 11])=topology(ind1,[5 7]);
   topology(ind2,[10 11])=topology(ind2,[6 8]);
   patch('faces',topology(:,[1 5 8]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle','none');
   hold on;
   patch('faces',topology(:,[5 2 6]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle','none');
   patch('faces',topology(:,[6 3 7]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle','none');
   patch('faces',topology(:,[7 4 8]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle','none');
   patch('faces',topology(:,[7 8 10]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle',':');
   patch('faces',topology(:,[5 6 11]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle',':');
   patch('faces',topology(:,1:4),'vertices',xyz,'FaceColor','none');
   hold off;
else
   patch('faces',topology(:,1:nknel),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp');
end
