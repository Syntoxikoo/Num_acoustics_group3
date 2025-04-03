function [nodesb,elementsb,segments]=meshcheck(nodes,elements,varargin)

% [nodesb,elementsb,segments]=meshcheck(nodes,elements,faulty,plotgeom,elsizes)
%
% Checks a 3D geometry produced by a meshing program. Finds body numbers
% and checks normal vectors. All normal vectors are set to point outwards.
% Checks if more than two elements are connected by an edge and if there
% are open edges, and signals it.
%
% Each body in the geometry is drawn in a separate figure, indicating the
% open edges in red and the illegally connected edges in blue.
%
% WARNINGS:
%  -"meshcheck" is a pre-check to detect the most common problems related
%  with BEM meshes. It can only fix normal vector directions and remove
%  unused nodes. Other errors, like badly connected or unconnected elements
%  are found, but not fixed.
%  -There are other possible mesh errors that "meshcheck" cannot find. Some
%  of these problems can provoke "meshcheck" to crash and others may run
%  unnoticed. For example, if two bodies are topologically independent but
%  physically share the same space, "meshcheck" will not give any error.
%  -The user will have to correct the mesh problems using a meshing program.
%
% Input variables:
%    -nodes:     node x, y and z coordinates. One node per row.
%                Do not include body numbers.
%    -elements:  numbers of the nodes in each element, ordered.
%                One element per row. Do not include body numbers.
%                It is possible to include elements with different 
%                number of nodes. The remaining of the shorter rows will
%                not be read and can be set, i.e., to NaN.
%    -faulty:    If this flag is set to 1, "meshcheck" will output the
%                nodes and elements matrices also for faulty meshes. if
%                this parameter is not included or is equal to 0, the 
%                matrices are set to "NaN" in case of errors.
%    -plotgeom:  Values: 0, no plot; 1, one plot per body (default); 2, the
%                plots are drawn element by element. The last option is
%                slow for large geometries, but it can be run with portions
%                of the elements matrix, in order to see details.
%    -elsizes:   a column vector containing the number of nodes in each
%                element. Its length is the number of elements. Can be
%                omitted if all elements have the same number of nodes.
%
% Output variables:
%    -nodesb:    node x, y and z coordinates and body number.
%                One node per row.
%    -elementsb: numbers of the nodes in each element, ordered,
%                plus body number. One element per row.
%    -segments:  matrix with one row per segment. A segment are two
%                connected nodes in an element. The columns are: node 1,
%                node 2, element 1, element 2, body number. If the segment
%                is connected to more than one element, column 4 is set to
%                Inf. The segment will then exist in more than one body.
%                If the segment is open (one element), column 4 is NaN.
%
% The element chosen to decide the overall direction is painted in a 
% darker shade: it is the farthest from the center of mass.
%
% The user must change manually the signs of every body number
% in the node list (nodesb) and in the element list (topologyb) to
% indicate the BEM which domain to consider. The convention is:
%                    +: exterior domain
%                    -: interior domain
% By default all normal vectors are set to point to the exterior domain (+).

% Vicente Cutanda Henríquez, 5-2011


if size(nodes,2)>3
   error('Error: do not include body numbers in the input arrays - Checking aborted');
end

[N,ncol]=size(elements);
elements=[elements ones(N,1)*NaN];
 
M=size(nodes,1);
nodes=[nodes ones(M,1)*NaN];

segments=[];

% Plot color definition
%Cdata=ones(N,1)*[0.7 .92 0.75]; % color data, B&K green
%Cdata=ones(N,1)*[209 85 34]/255; % color data, D-A-S orange
%Cdata=ones(N,1)*[204 204 204]/255; % color data, D-A-S light grey
Cdata=ones(N,1)*[233 229 218]/255; % color data, D-A-S beige
%Cdata=ones(N,1)*[68 68 136]/255; % color data, D-A-S dark blue
ColorOpen=[1 0 0];ColorIllegal=[0 0 1];

% Fill input variables
if nargin==2
    faulty=0;
elseif nargin>2
    faulty=varargin{1};
end

if nargin<4
    plotgeom=1;
else
    plotgeom=varargin{2};
end

if nargin<5
    elsizes=ones(N,1)*ncol; % All elements are of the same number of nodes, no need to provide "elsizes"
else
    elsizes=varargin{3};
end


% Check geometry body by body
body=0;
while any(isnan(elements(:,end))) % body loop until all elements have been checked
   body=body+1;
   
   elnr=find(isnan(elements(:,end)));elnr=elnr(1);
   chg_sgn=0;
   % Call "elstep" function. It sweeps over the body surface, checking
   % the connections between elements, changing normal vectors if needed, 
   % assigning body numbers and reporting open segments and cases where
   % more than two elements are connected with a segment.
   if plotgeom==2 % pass color data for the plots if a element by element plot is required
       [segments,nodes,elements]=elstep(segments,nodes,elements,elnr,body,elsizes,chg_sgn,Cdata,ColorOpen,ColorIllegal);
   else
       [segments,nodes,elements]=elstep(segments,nodes,elements,elnr,body,elsizes,chg_sgn);
   end
end   


% Remove the unused nodes and adjust the matrices.
if any(isnan(nodes(:,end)))
    shifts=zeros(M,1);
    for jj=1:M % Find how much to shift every node
        if isnan(nodes(jj,end))
            shifts(jj)=NaN;
            shifts(jj+1:end)=shifts(jj+1:end)+1;
        end
    end
    for ff=1:N % Shift every node reference in "elements"
        for cc=1:elsizes(ff)
            elements(ff,cc)=elements(ff,cc)-shifts(elements(ff,cc));
        end
    end
    for ff=1:size(segments,1) % Shift every node reference in "segments"
        segments(ff,1)=segments(ff,1)-shifts(segments(ff,1));
        segments(ff,2)=segments(ff,2)-shifts(segments(ff,2));
    end
    nodes=nodes(~isnan(shifts),:); % Remove the unused nodes from the list
    disp([num2str(sum(isnan(shifts))) ' unused nodes were removed from the node list'])
end


for bb=1:body
   tx=find(elements(:,end)==bb);
   Nb=length(tx);
   nodestmp=nodes(nodes(:,end)==bb,1:3);

   %%%%%%%% Check if all normal vectors point outwards, and change if not
   
   cent_b=mean(nodestmp); % mass center of the body
   % find the farthest element
   dist=0;
   for jj=1:Nb
       nnod=elsizes(tx(jj));
       dist_temp=mean(sqrt(sum((nodes(elements(tx(jj),1:nnod),1:3)-repmat(cent_b,nnod,1)).^2,2))); % mean distance to the nodes
       if dist_temp > dist
           dist=dist_temp;
           far_ele=tx(jj);
       end
   end
   Cdata(far_ele,:)=Cdata(far_ele,:)*0.6; %indicate the farthest element by a color change
   
   farnod=nodes(elements(far_ele,1:3),1:3); % coordinates of the three first nodes of the farthest element
   % Take the cross product of the two first segments of the far element to
   % get a normal vector at its second node. Then examine the sign of the
   % dot product of the result with the vector from the second node to the
   % mass center. If the sign is negative, the normal vector points to the
   % interior of the body.
   if dot(cross(farnod(3,:)-farnod(2,:),farnod(1,:)-farnod(2,:)),cent_b-farnod(2,:))>0
       % change all normal vectors in the body if they point inwards
       for jj=1:Nb
           nnod=elsizes(tx(jj));
           elements(tx(jj),1:nnod)=elements(tx(jj),nnod:-1:1);
       end
   end
   
   %%%%%%%% Plot geometry
   
   if plotgeom==1
       figure; % one figure for each body
       for jj=1:Nb % draw elements in the current body one by one
           nnod=elsizes(tx(jj));
           patch('faces',1:nnod,'vertices',nodes(elements(tx(jj),1:nnod),1:3),'FaceVertexCData',ones(nnod,1)*Cdata(tx(jj),:),'FaceColor','flat');
           hold on;
       end
%        % plot normal vectors on the second node of every element
%        normvects=cross(nodes(elements(tx,3),1:3)-nodes(elements(tx,2),1:3),...
%            nodes(elements(tx,1),1:3)-nodes(elements(tx,2),1:3));
%        quiver3(nodes(elements(tx,2),1),nodes(elements(tx,2),2),nodes(elements(tx,2),3),...
%            normvects(:,1),normvects(:,2),normvects(:,3));
       view(-135,45);grid
       rotate3d on;
       title(['Body nr. ' num2str(bb)]);
       xlabel('x');ylabel('y');zlabel('z');
       
       segbody=segments(segments(:,end)==bb,:);
       for jj=1:size(segbody,1)
           % Draw red lines indicating open segments
           if isnan(segbody(jj,4))
               h=line(nodes(segbody(jj,1:2),1),nodes(segbody(jj,1:2),2),nodes(segbody(jj,1:2),3));
               set(h,'LineWidth',1.5,'Color',ColorOpen);
           end
           % Draw blue lines indicating illegal connections with other bodies
           if isinf(segbody(jj,4))
               h=line(nodes(segbody(jj,1:2),1),nodes(segbody(jj,1:2),2),nodes(segbody(jj,1:2),3));
               set(h,'LineWidth',1.5,'Color',ColorIllegal);
           end
       end
       hold off;
   end
   
end

% set output matrices to NaN if there are errors
if any(isnan(segments(:,4))) || any(isinf(segments(:,4)))
    disp('Errors were detected in the mesh. It is not suitable for BEM calculations.')
    if faulty
        nodesb=nodes;
        elementsb=elements;
    else
        nodesb=NaN;
        elementsb=NaN;
    end
else
    disp('No topological errors were detected in the mesh. All normal vectors are set to point outwards.')
    nodesb=nodes;
    elementsb=elements;
end

end


function [segments,nodes,elements]=elstep(segments,nodes,elements,elnr,body,elsizes,chg_sgn,varargin)

% Analyzes an element's segments, reorders if needed. Finds the
% elements attached, evaluates if the connection is good, and runs
% again for the new elements not yet analyzed. It does so until a body has
% been analyzed.
%
% All normal vectors are set to point in the same direction.
%
% Input variables:
%    -segments:  matrix with one row per segment. A segment are two
%                connected nodes in an element. The columns are: node 1,
%                node 2, element 1, element 2, body number. If the segment
%                is connected to more than one element, column 4 is set to
%                Inf. The segment will then exist in more than one body.
%                If the segment is open (one element), column 4 is NaN.
%    -nodes:     node x, y and z coordinates. One node per row.
%                Body numbers will be added.
%    -elements:  numbers of the nodes in each element, ordered.
%                One element per row. Body numbers are included in an extra
%                column. It is possible to include elements with different 
%                number of nodes. The remaining of the shorter rows will
%                not be read and can be set, i.e., to NaN.
%    -elnr:      number of the element to be analyzed
%    -body:      body number inherited from the first call
%    -elsizes:   a column vector containing the number of nodes in each
%                element. Its length is the number of elements. 
%                A minus sign indicates that the normal vector sign was
%                changed.
%    -chg_sgn:   if 1, the order of the nodes in the element is changed
%                (sign of the normal vector).
%
% Output variables: segments, nodes and elements, modified.

% Vicente Cutanda Henríquez, 5-2011

ele_stack=elnr;

if nargin>7
    Cdata=varargin{1};ColorOpen=varargin{2};ColorIllegal=varargin{3};
    % Initialize element by element figure
    figure; % one figure for each body
    hold on;
    view(-135,45);grid
    rotate3d on;
    title(['Body nr. ' num2str(body)]);
    xlabel('x');ylabel('y');zlabel('z');
end

while ~isempty(ele_stack)
    
    % take one element and remove it from the stack
    elnr=ele_stack(1); ele_stack=ele_stack(2:end);
    
    elements(elnr,end)=body;
    nnod=abs(elsizes(elnr));
    elek=elements(elnr,1:nnod);

    if nargin>7
        patch('faces',1:nnod,'vertices',nodes(elements(elnr,1:nnod),1:3),'FaceVertexCData',ones(nnod,1)*Cdata(elnr,:),'FaceColor','flat'); drawnow
%         % plot normal vectors on all nodes of every element
%         for nn=1:nnod
%             normvects=cross(nodes(elements(elnr,mod(nn,nnod)+1),1:3)-nodes(elements(elnr,nn),1:3),...
%                 nodes(elements(elnr,mod(nn-2,nnod)+1),1:3)-nodes(elements(elnr,nn),1:3));
%             quiver3(nodes(elements(elnr,nn),1),nodes(elements(elnr,nn),2),nodes(elements(elnr,nn),3),...
%                 normvects(:,1),normvects(:,2),normvects(:,3));
%         end
    end
    
    for ee=1:nnod % examine all segments in the element
        n1=elek(ee);n2=elek(mod(ee,nnod)+1); % nodes in the current segment
        
        % Determine if the segment is already listed in the current body
        if ~isempty(segments)
            inlist=find(((segments(:,1)==n1 & segments(:,2)==n2) | (segments(:,2)==n1 & segments(:,1)==n2) | ...
                (segments(:,1)==n2 & segments(:,2)==n1)) & segments(:,5)==body);
        else
            inlist=[];
        end
        
        if length(inlist)>1
            error(['Calculation error,  element: ' num2str(elnr) ', nodes: ' num2str(n1) ', ' num2str(n2)]);
        elseif isempty(inlist)
            %       elnr_new=[];segnr_new=[];
            [ff,cc]=find(elements(:,1:end-1)==n1); % find new elements with the same segment
            ffel=find(ff==elnr);
            
            if length(ffel)==1 % exclude the element under examination
                ff=ff([1:ffel-1 ffel+1:end]);cc=cc([1:ffel-1 ffel+1:end]);
            else
                error(['Meshcheck error. Possible causes: element with repeated nodes, element: ' ...
                    num2str(elnr) ', nodes: ' num2str(n1) ', ' num2str(n2)]);
            end
            
            ne=0;ff2=[];cc2=[];
            for ffi=1:length(ff) % Find other elements sharing the segment.
                if elements(ff(ffi),mod(cc(ffi),nnod)+1)==n2
                    ff2=[ff2 ff(ffi)];cc2=[cc2 cc(ffi)];
                    chg_sgn=1; % new element found, with normal vector in the opposite direction
                    ne=ne+1;
                elseif elements(ff(ffi),mod(cc(ffi)-2,nnod)+1)==n2
                    ff2=[ff2 ff(ffi)];cc2=[cc2 cc(ffi)];
                    chg_sgn=0; % new element found, with normal vector in the same direction
                    ne=ne+1;
                end % if only one node is shared ("else" cases), do nothing
            end
            
            if ne>1
                disp(['More than two elements are connected to a common segment! Element: ' ...
                    num2str(elnr) ', nodes: ' num2str(n1) ', ' num2str(n2)]);
                segments=[segments; n1 n2 elnr Inf body];
                nodes(n1,end)=body; nodes(n2,end)=body;
                if nargin>7
                    h=line(nodes([n1 n2],1),nodes([n1 n2],2),nodes([n1 n2],3));
                    set(h,'LineWidth',1.5,'Color',ColorIllegal);
                end
            elseif ne==0
                disp(['There is an open edge! Element: ' num2str(elnr) ', nodes: ' num2str(n1) ', ' num2str(n2)])
                segments=[segments; n1 n2 elnr NaN body];
                nodes(n1,end)=body; nodes(n2,end)=body;
                if nargin>7
                    h=line(nodes([n1 n2],1),nodes([n1 n2],2),nodes([n1 n2],3));
                    set(h,'LineWidth',1.5,'Color',ColorOpen);
                end
            else % a new segment has been found and the neighboring element is added to the stack for later analysis
                segments=[segments; n1 n2 elnr ff2 body];
                nodes(n1,end)=body; nodes(n2,end)=body;
                if any(ele_stack==ff2)
                    if chg_sgn
                        error(['Non-circular node ordering in element: ' num2str(elnr) ', nodes: ' num2str(n1) ', ' num2str(n2)]);
                    end
                else
                    ele_stack=[ele_stack ; ff2];
                    if chg_sgn % invert order of the nodes if needed
                        elements(ff2,1:nnod)=elements(ff2,nnod:-1:1);
                    end
                end
            end
        end  % do nothing if length(inlist)==1, when the segment is already listed
        
    end
end
end

