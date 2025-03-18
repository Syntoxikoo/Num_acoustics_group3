function [NodeRhoZ,ElemNodeNum,M,N,NumBodies,DataBody]=readnods(File,see)

%  
%  [NodeRhoZ,ElemNodeNum,M,N,NumBodies,DataBody]=readnods(File,see)
%  
%  Reads a text file that contains node coordinates, as it is written by
%  'nodegen', and places data into variables:
%  
%    -NodeRhoZ : Node coordinates along the generator. Matrix with M rows,
%      first column is the node rho coodinate, second column is its z
%      coordinate and third column is the number of the body it belongs to,
%      with a minus sign if the interior domain is specified.
%    -ElemNodeNum : Contains the node numbers for the 3 nodes within each
%      element. Four columns and N rows. The fourth column has the
%      corresponding body number, with a minus sign if the interior domain
%      is specified. 
%    -M : Total number of nodes. Single integer.
%    -N : total number of elements. Single integer.
%    -NumBodies : Number of separate bodies. Single integer.
%    -DataBody : Matrix with NumBodies rows, each row contains one body data.
%      First column is number of elements, second is the number of nodes,
%      third is the starting element number, fourth is the last element number,
%      fifth is the starting node number, sixth is the last node number and
%      seventh is 1 if the domain is exterior or -1 if it is interior.
%  
%  If 'see' is set to 'y' or 'Y', a plot is made with the data. Each body is
%  plotted blue if the domain is exterior to it, or red if it is interior.
%  
%  
%  Notes:
%    -Diaphragm data is not read in this version.
%    -There may be more data than necessary for the calculations.
%    -The input file must have must have minus sign on the number of elements
%     for the bodies to which the domain is interior.

%  Vicente Cutanda 1998
%  VC: Modified to add signs for interior/exterior domains (2.2000)
%  

nodepos=fopen(File);
AA=fscanf(nodepos,'%u %u',[2]);
NumBodies=AA(1,1);
M=0;N=0;
for i=1:NumBodies
  AA=fscanf(nodepos,'%g %g',[2]);
  DataBody(i,1)=abs(AA(1,1));
  interior=sign(AA(1,1));
  for j=(N+1):(N+DataBody(i,1))
    ElemNodeNum(j,1)=2*j+i-2;
    ElemNodeNum(j,2)=2*j+i-1;
    ElemNodeNum(j,3)=2*j+i;
    ElemNodeNum(j,4)=i*interior;
  end

  DataBody(i,3)=N+1;
  N=N+DataBody(i,1);
  DataBody(i,4)=N;

  DataBody(i,2)=2*DataBody(i,1)+1;

  DataBody(i,5)=M+1;
  M=M+DataBody(i,2);
  DataBody(i,6)=M;
  DataBody(i,7)=interior;
end

for i=1:M
  AA=fscanf(nodepos,'%g %g',[2]);
  AA=AA';
  NodeRhoZ(i,1:2)=AA(1,1:2);
  for t=1:NumBodies
    if i>=DataBody(t,5) & i<=DataBody(t,6)
      NodeRhoZ(i,3)=t*DataBody(t,7);
    end
  end
end
fclose(nodepos);

if see(1)=='y' | see(1)=='Y'
   for bb=1:NumBodies
      if DataBody(bb,7)==1
         plot(NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),1),1)...
             ,NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),1),2),'ob');
         hold on;
         plot(NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),2),1)...
             ,NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),2),2),'+b');
         plot(NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),3),1)...
             ,NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),3),2),'ob');
      else
         plot(NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),1),1)...
             ,NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),1),2),'or');
         hold on;
         plot(NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),2),1)...
             ,NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),2),2),'+r');
         plot(NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),3),1)...
             ,NodeRhoZ(ElemNodeNum(DataBody(bb,3):DataBody(bb,4),3),2),'or');
      end
   end
    
   hold off;
   title(['Nodes = ' int2str(M) '  Elements = ' int2str(N) '  Bodies = ' int2str(NumBodies)]);
   xlabel('rho');
   ylabel('z');
   Xmax=max(NodeRhoZ(:,1));
   Ymax=max(abs(NodeRhoZ(:,2)));
   Xmax=Xmax+Xmax*0.3;
   Ymax=Ymax+Ymax*0.5;
   axis([0 Xmax -Ymax Ymax]);
end
