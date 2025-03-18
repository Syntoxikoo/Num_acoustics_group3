function [nodes,elementsTRI,elementsQUAD]=readgeomGMSH(filename)

% [nodes,elementsTRI,elementsQUAD]=readgeomGMSH(filename);
%
% reads nodes (x,y,z) and elements into matrices
% Version able to read from a GMSH format file

[fid,message] = fopen(filename);

fin=0;
while feof(fid) == 0 & fin==0
   line=fgetl(fid);
   if length(line)==6
       if line=='$Nodes'
           fin=1;
       end
   end
end
line=fgetl(fid);
M=sscanf(line,'%d');

for jj=1:M
    line=fgetl(fid);
    goodline=sscanf(line,'%d%f%f%f');
    nodes(jj,:)=[goodline(2) goodline(3) goodline(4)];
end

fin=0;
while feof(fid) == 0 & fin==0
   line=fgetl(fid);
   if length(line)==9
       if line=='$Elements'
           fin=1;
       end
   end
end
line=fgetl(fid);
N=sscanf(line,'%d');
elementsTRI=[];elementsQUAD=[];
fin=0;
while feof(fid) == 0 & fin==0
    line=fgetl(fid);
    if length(line)==12 & line=='$EndElements'
        fin=1;
    else
        goodline=sscanf(line,'%d');
        if goodline(2)==2 % 3-node triangle
            elementsTRI=[elementsTRI ; goodline(end-2:end)'];
        end
        if goodline(2)==9 % 6-node triangle
            elementsTRI=[elementsTRI ; goodline(end-5:end)'];
        end
        if goodline(2)==3 % 4-node quadrangle
            elementsQUAD=[elementsQUAD ; goodline(end-3:end)'];
        end
        if goodline(2)==21 % 10-node cubic triangle
            elementsTRI=[elementsTRI ; goodline(end-9:end)'];
        end
    end
end

fclose(fid);
