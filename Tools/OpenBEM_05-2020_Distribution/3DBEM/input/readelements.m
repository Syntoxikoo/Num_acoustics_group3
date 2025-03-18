function topology=readelements(filename)

% topology=readelements(filename)
% reads corner nodes global numbers into Nx4 matrix
% N number of elements

[fid,message] = fopen(filename);
elemnum=1;
topology=[];
while feof(fid) == 0
   line=fgetl(fid);
   linestart=sscanf(line,'%s',1);
   if ~isempty(linestart) & ~isletter(linestart)
      if eval(linestart)==elemnum
         goodline=sscanf(line,'%d');
         l=length(goodline);
         topology=[topology; goodline(7:end)'];
         elemnum=elemnum+1;
      end
      
   end
   
end
fclose(fid);
%rearrange all normals since ANSYS favours opposite normals than I
all=1:size(topology,1);
topology=rearrange(topology,all);