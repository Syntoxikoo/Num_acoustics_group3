function [BCtopo,BCnodeA,BCnodeB]=bound(xyzb,topology,BCelem,see)
% [BCtopo,BCnodeA,BCnodeB]=bound(xyzb,topology,BCelem,see);
%
% Creates matrices to handle expanded columns in the B matrix
% due to element-oriented boundary conditions. The surface is
% divided into numbered boundary areas with continuous boundary
% conditions.
%
% Input variables:
%    -xyzb:      node x, y and z coordinates. One node per row.
%                Include body numbers.
%    -topology:  numbers of the nodes in each element, ordered.
%                One element per row. Include body numbers.
%    -BCelem:    column vector with N (number of elements) values.
%                Each value is the boundary area number the
%                element belongs to (1,2,3...).
%    -see:       if "y", plots the geometry with boundary areas.
%                The areas are darker for higher area numbers and
%                the nodes on area boundaries are drawn in red.
%
% Output variables:
%    -BCtopo:    same structure as topology, but it contains
%                indexes (rows) in BCnodeB instead of node numbers.
%    -BCnodeA:   matrix with M (number of nodes) rows. It has as
%                many columns as boundary areas and the elements
%                are 1 if the node (row number) belongs to the
%                area (column number) or 0 if it does not.
%    -BCnodeB:   matrix with two columns: the first are node
%                numbers, and the second are boundary area
%                numbers they belong to. If a node belongs to
%                several areas, there is one row per area.
%                Rows correspond to columns in the B matrix.

% Vicente Cutanda 09.2000


[N,nknel] = size(topology);
topology=abs(real(topology));
nknel=nknel-1;
[M,cols]=size(xyzb);
if cols<4
   error('Error: The input arrays must include body numbers - Calculation aborted');
end
xyzb(:,4)=abs(xyzb(:,4));

BCnodeA=zeros(M,max(BCelem));
BCnodeB=[];
for bb=1:M
   [rr,cc]=find(topology(:,1:nknel)==bb); % Find elements where node bb belongs
   [tmp,ind]=sort(BCelem(rr)); % Get  and sort zones these elements belong to, and link to the elements (ind). Zones may be repeated.
   ind2=[1 ; find(diff(tmp)~=0)+1]; % Select only indexes where these is a zone change
   tmp2=tmp(ind2); % Zones again, but only the non-repeated.
   ll=size(BCnodeB,1); % Count B columns so far
   BCnodeB=[BCnodeB ; bb*ones(length(tmp2),1) tmp2]; 
   BCnodeA(bb,tmp2)=1;
   
   rr=rr(ind);cc=cc(ind); % Sort the element positions in topology as growing zone number
   tmp3=abs(sign([1;diff(tmp)])); 
   aa=0;
   for tt=1:length(tmp3)
      aa=aa+tmp3(tt);
      BCtopo(rr(tt),cc(tt))=aa+ll;
   end
end

% draw plot if requested
if see=='y'
   NumBodies = max(xyzb(:,4));
   for bb=1:NumBodies
      tx=find(topology(:,end)==bb);
      Nb=length(tx);
      nx=find(xyzb(:,end)==bb);
      if nknel==8
         colstopo=[1 5 2 6 3 7 4 8];
      elseif nknel==6
         colstopo=[1 4 2 5 3 6];
      else
         colstopo=1:nknel;
      end
      topologytmp=topology(tx,colstopo);
      nodestmp=xyzb(nx,1:3);
      % color data, shades of B&K green 
      Cdata=(max(BCelem(tx))+2-BCelem(tx))./max(BCelem(tx)+1)*[0.7 .92 0.75];

      figure; % one figure for each body, in case they hide one another
      view(-135,45);
%      axis equal;
      rotate3d on;
      hold on;
      title(['Body nr. ' num2str(bb)]);
      xlabel('x');ylabel('y');zlabel('z');
      patch('faces',topologytmp,'vertices',nodestmp,'FaceVertexCData',Cdata,'FaceColor','flat');
   end
   
   if nknel==8 % plot nodes 5 to 8
      temp=topology(:,5:8);
      plot3(xyzb(temp(1:end),1),xyzb(temp(1:end),2),xyzb(temp(1:end),3),'k.');
   end
   if nknel==6 % plot nodes 4 to 6
      temp=topology(:,4:6);
      plot3(xyzb(temp(1:end),1),xyzb(temp(1:end),2),xyzb(temp(1:end),3),'k.');
   end
   [rrr,tmp]=find(sum(BCnodeA')'>1); 
   plot3(xyzb(rrr,1),xyzb(rrr,2),xyzb(rrr,3),'r.'); % plot boundary nodes in red
   
end
