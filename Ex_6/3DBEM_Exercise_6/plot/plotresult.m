function []=plotresult(xyzb,topology,result,varargin)
% []=plotresult(xyzb,topology,result,mode)
% Plots result as color on a patch surface
%
% xyzb : nodal coordinates
% topology: relation between elements and nodes
% result : vector to be plotted
% mode : 1, modulus; 2, phase; 3, both. Default: real input

[N,nknel] = size(topology);
[M,cols]=size(xyzb);
if cols~=3
   xyz=xyzb(:,1:3);
   nknel=nknel-1;
else
   xyz=xyzb;
end


figure;

if nargin>3
   mode=varargin{1};
   if mode==1
      drawone(xyz,topology,abs(result),nknel)
   elseif mode==2
      drawone(xyz,topology,angle(result)*180/pi,nknel)
   elseif mode==3
      subplot(1,2,1)
      drawone(xyz,topology,abs(result),nknel)
      colorbar vert;
      xlabel('x');ylabel('y');zlabel('z');
      subplot(1,2,2)
      drawone(xyz,topology,angle(result)*180/pi,nknel)
   end
else
   drawone(xyz,topology,result,nknel)
end

colorbar vert;
xlabel('x');ylabel('y');zlabel('z');
%axis equal;
grid
rotate3d on;


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
   patch('faces',topology(:,[7 8 10]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle','none');
   patch('faces',topology(:,[5 6 11]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp','LineStyle','none');
   patch('faces',topology(:,1:4),'vertices',xyz,'FaceColor','none');
   hold off;
else
   patch('faces',topology(:,1:nknel),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp');
end
