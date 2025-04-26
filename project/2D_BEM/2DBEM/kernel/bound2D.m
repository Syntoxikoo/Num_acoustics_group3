function [BCtopo,BCnodeA,BCnodeB]=bound2D(xyb,topology,BCelem,see)
% [BCtopo,BCnodeA,BCnodeB]=bound(xyb,topology,BCelem,see);
%
% Creates matrices to handle expanded columns in the B matrix
% due to element-oriented boundary conditions. The surface is
% divided into numbered boundary areas with continuous boundary
% conditions. Version in two dimensions.
%
% Input variables:
%    -xyb:       node x and y coordinates. One node per row.
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

% Vicente Cutanda 09.2000, original in 3DBEM
% Version in two dimensions, VCH 4-2017


[N,nknel] = size(topology);
topology=abs(real(topology));
nknel=nknel-1;
[M,cols]=size(xyb);
if cols<3
   error('Error: The input arrays must include body numbers - Calculation aborted');
end
xyb(:,3)=abs(xyb(:,3));

BCnodeA=zeros(M,max(BCelem));
BCnodeB=[];
for bb=1:M
   [rr,cc]=find(topology(:,1:nknel)==bb); % Find elements where node bb belongs
   [tmp,ind]=sort(BCelem(rr)); % Get and sort zones these elements belong to, and link to the elements (ind). Zones may be repeated.
   ind2=[1 ; find(diff(tmp)~=0)+1]; % Select only indexes where there is a zone change
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
    
    figure;
    for tt=1:size(topology,1) % Plot segments alternating colors
        if nknel==3
            switch BCelem(tt)
                case 1, colel='r'; case 2, colel='b'; case 3, colel='g'; case 4, colel='y'; case 5, colel='m';
            end
            plot(xyb(topology(tt,1:3),1),xyb(topology(tt,1:3),2),[colel ':']);
            hold on;
            plot(xyb(topology(tt,[1 3]),1),xyb(topology(tt,[1 3]),2),[colel 'o']);
            plot(xyb(topology(tt,2),1),xyb(topology(tt,2),2),[colel '+']);
        elseif nknel==2
            switch BCelem(tt)
                case 1, colel='r'; case 2, colel='b'; case 3, colel='g'; case 4, colel='y'; case 5, colel='m';
            end
            plot(xyb(topology(tt,1:2),1),xyb(topology(tt,1:2),2),[colel ':o']);
            hold on;
%            plot(xyb(topology(tt,1:2),1),xyb(topology(tt,1:2),2),[colel 'o']'ko');
        end
            
    end
    plot(xyb(sum(BCnodeA,2)>1,1),xyb(sum(BCnodeA,2)>1,2),'r*')
    
    hold off;
    grid;
    axis equal;
    title(['Nodes = ' num2str(size(xyb,1)) '  Elements = ' num2str(size(topology,1)) ...
           ' Split nodes = ' num2str(size(BCnodeB,1)-size(xyb,1)) '  Boundary regions = ' num2str(max(BCnodeB(:,2)))]);
    xlabel('x');
    ylabel('y');
    
end
