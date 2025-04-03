function nodes=readnodes(filename)

% nodes=readnodes(filename)
% reads nodes (x,y,z) into Mx3 matrix
% M number of nodes

[fid,message] = fopen(filename);
nodenum=1;
nodes=[];
while feof(fid) == 0
   line=fgetl(fid);
   linestart=sscanf(line,'%s',1);
   if ~isempty(linestart) & ~isletter(linestart)
      if eval(linestart)==nodenum
         goodline=sscanf(line,'%d%f%f%f');
         nodes=[nodes; goodline(2) goodline(3) goodline(4)];
         nodenum=nodenum+1;
      end
      
   end
   
end
%figure;
%plot3(nodes(:,1),nodes(:,2),nodes(:,3),'*')
fclose(fid);