function []=plotresult(xyzb,topology,result)
% []=plotresult(xyzb,topology,result)
% Plots result as color on a patch surface
%
% xyzb : nodal coordinates
% topology: relation between elements and nodes
% result : real vector to be plotted

[N,nknel] = size(topology);
[M,cols]=size(xyzb);
if cols~=3
   xyz=xyzb(:,1:3);
   nknel=nknel-1;
else
   xyz=xyzb;
end

figure;
if nknel==8
   for tt=1:N
      %xyz(M+tt,1:3)=[mean(xyz(topology(tt,1:nknel),1)) mean(xyz(topology(tt,1:nknel),2)) ...
      %      mean(xyz(topology(tt,1:nknel),3))];
      [psi, xq, yq, zq, nx, ny, nz]=elemshape(xyzb(topology(tt,1:nknel),:),[0 0 0 0]);
      xyz(M+tt,1:3)=[xq yq zq];
      topology(tt,9)=M+tt;
      result(M+tt)=mean(result(topology(tt,1:nknel)));
   end
   patch('faces',topology(:,[8 4 7 9]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp');
   hold on;
   patch('faces',topology(:,[9 7 3 6]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp');
   patch('faces',topology(:,[5 9 6 2]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp');
   patch('faces',topology(:,[1 8 9 5]),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp');
   plot3(xyz(M+1:end,1),xyz(M+1:end,2),xyz(M+1:end,3),'k.');
else
   patch('faces',topology(:,1:nknel),'vertices',xyz,'FaceVertexCData',result,'FaceColor','interp');
end


colorbar vert;
xlabel('x');ylabel('y');zlabel('z');
%axis equal;
rotate3d on;
