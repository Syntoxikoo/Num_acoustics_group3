 function [nodesb,topologyb,toposhrinkb,tim,segmopen,segmextraele]=stlbodyfind(nodes,topology)

% [nodesb,topologyb,toposhrinkb,tim,segmopen]=stlbodyfind(nodes,topology);
% modified version of Vicentes bodyfind
% Changes: 
%          uses elemshapestl (triangular elements), which is now 'built-in'
% revised: 11. dec 2006 pmj


% Finds body numbers and checks normal vectors.
% All normal vectors are set to point outwards.
%
% Input variables:
%    -nodes:     node x, y and z coordinates. One node per row.
%                Do not include body numbers.
%    -topology:  numbers of the nodes in each element, ordered.
%                One element per row. Do not include body numbers.
%
% Output variables:
%    -nodesb:    node x, y and z coordinates and body number.
%                One node per row.
%    -topologyb: numbers of the nodes in each element, ordered,
%                plus body number. One element per row.
%    -toposhrinkb:numbers of the nodes in each element, ordered,
%                plus body number. Only the elements that had a
%                duplicated node are included, with that node
%                supressed. One element per row.
%    -tim:       calculation time.
%    -segmopen:  list of segments that are open edges. First column
%                is the element it belongs to, second column is
%                the first node of the segment within the element,
%                and third column is the body number.
%    -segmextraele: list of segments belonging to more than two elements.
%                First column is the element it belongs to, second column
%                is the first node of the segment within the element,
%                and third column is the body number.(VCH addition 08/2008)
%
% Each body is plotted in a separate figure, a different shade
% of green indicating which elements had an opposite direction.
% The element chosen to decide the overall direction is painted
% in darker green: it is the farthest from the center of mass.
% The open edges, if any, are indicated in blue.
%
% The user must change manually the signs for every body number
% in the node list (nodesb) and in the element list (topologyb) to
% indicate the BEM which domain to consider. The convention is:
%                    +: exterior domain
%                    -: interior domain

% Vicente Cutanda 08.2000
% Peter Juhl 08.2001 revised aug 2008


tic

if size(nodes,2)>3
   error('Error: do not include body numbers in the input arrays - Checking aborted');
end

% Removed from original bodyfind
[N,nknel]=size(topology);
isquad=0;
% if nknel==8
%    nknel=4;
%    isquad=1;
% else
%    isquad=0;
% end

topology=[topology zeros(N,1)];

segmopen=[];
elemcorr=[];
segmextraele=[];
body=0;
while ~isempty(find(topology(:,end)==0)) % body loop
    body=body+1
    elez=find(topology(:,end)==0);
    eln=1;
    segms=[elez(ceil(eln/4)) mod((eln-1),nknel)+1];  % start a stack list of segments to test
    flag=1; % first segment is different

    while ~isempty(segms) % segment loop
        randseg=1; %ceil(rand*size(segms,1));   %it can be random
        segmin=segms(randseg,:); % choose one segment in the list
        segms(randseg,:)=[]; % remove it from the list

        % changed from quad to isquad in last argument
        [segmout,topology]=elstep(topology,segmin,body,nknel,N,isquad);

        if size(segmout,1)>1

            % START: newly added code to cope with more than two adjoining elements (VCH 08/2008)
            disp(['There are more than two elements with a a segment in common! - Element nr.: ' ...
                        num2str(segmin(1)) ' Segment pos.: ' num2str(segmin(2))]);
            segmextraele=[segmextraele ; segmin body];
            for jj=1:size(segmout,1)
                if segmout(jj,1)<0
                    segmout(jj,:)=abs(segmout(jj,:));
                    disp(['Opposite element direction found and corrected - Element nr.: ' ...
                          num2str(segmout(jj,1))]);
                    elemcorr=[elemcorr ; segmout(jj,1) body];
                end

                if topology(segmout(jj,1),end)==0 % unchecked element found
                    segms=[segms ; [ones(nknel+flag-1,1)*segmout(jj,1) ...
                        (mod((segmout(jj,2)-flag):(nknel+segmout(jj,2)-2),nknel)+1)']];

                    topology(segmout(jj,1),end)=body;
                    nodes(topology(segmout(jj,1),1:nknel*(1+isquad)),4)=body;
                elseif topology(segmout(jj,1),end)==body % checked element found
                    pos=find(sum((abs(segms-ones(size(segms,1),1)*segmout(jj,:)))')==0);
                    if length(pos)==1 % the adjoining segment should be in the stack list: remove it
                        segms(pos,:)=[];
                    else
                        error('Program error: segment count failure - Checking aborted');
                    end
                else % the new element belongs to another body!
                    error('Error: two bodies are connected - Checking aborted');
                end
            end
            flag=0;
            % END: newly added code to cope with more than two adjoining elements (VCH 08/2008)
            
        else
            if segmout(1)==0 % the segment is open, just add it to a list
                if flag==1 % just in case the first segment is open, get another instead
                    %disp('something wrong');
                    eln=eln+1;
                    segms=[elez(ceil(eln/4)) mod((eln-1),nknel)+1];
                else
                    disp(['There is an open edge! - Element nr.: ' ...
                        num2str(segmin(1)) ' Segment pos.: ' num2str(segmin(2))]);
                    segmopen=[segmopen ; segmin body];
                end
            elseif segmout(2)==0 % a node number is duplcated within an element, do nothing
                if flag==1 % just in case the first segment is duplicated, get another instead
                    eln=eln+1;
                    segms=[elez(ceil(eln/4)) mod((eln-1),nknel)+1];
                end
            else
                if segmout(1)<0
                    segmout=abs(segmout);
                    disp(['Opposite element direction found and corrected - Element nr.: ' ...
                          num2str(segmout(1))]);
                    elemcorr=[elemcorr ; segmout(1) body];
                end

                if topology(segmout(1),end)==0 % unchecked element found
                    segms=[segms ; [ones(nknel+flag-1,1)*segmout(1) ...
                        (mod((segmout(2)-flag):(nknel+segmout(2)-2),nknel)+1)']];

                    topology(segmout(1),end)=body;
                    nodes(topology(segmout(1),1:nknel*(1+isquad)),4)=body;
                elseif topology(segmout(1),end)==body % checked element found
                    pos=find(sum((abs(segms-ones(size(segms,1),1)*segmout))')==0);
                    if length(pos)==1 % the adjoining segment should be in the stack list: remove it
                        segms(pos,:)=[];
                    else
                        error('Program error: segment count failure - Checking aborted');
                    end
                else % the new element belongs to another body!
                    error('Error: two bodies are connected - Checking aborted');
                end
                flag=0;
            end
        end
    end
end


for bb=1:body
   tx=find(topology(:,end)==bb);
   Nb=length(tx);
   nx=find(nodes(:,end)==bb);
   topologytmp=topology(tx,1:nknel*(isquad+1));
   nodestmp=nodes(nx,1:3);
   
   
   % Cdata=[ones(Nb,1)]*[1 1 1]; %%%%%%%%%%
   Cdata=[ones(Nb,1)]*[0.7 .92 0.75]; % color data, B&K green
%   Cdata=[ones(Nb,1)]*[1 0 0]; % color data, spiderman version

   if ~isempty(elemcorr) %darken the elements that were corrected before
      corrbod=elemcorr(find(elemcorr(:,end)==bb),1);
      for ee=1:length(corrbod)
         %Cdata(find(tx==corrbod(ee)),:)=Cdata(find(tx==corrbod(ee)),:); %%%%%%%%%%
         Cdata(find(tx==corrbod(ee)),:)=Cdata(find(tx==corrbod(ee)),:)*0.9;
      end
   end
   
   for qq=1:length(nx)
      [ii,jj]=find(topology(tx,1:nknel)==nx(qq));
      for iii=1:length(ii)
         topologytmp(ii(iii),jj(iii))=qq;
      end
   end
   % NB should be modified for triangular elements pmj
   IP=[1/3 1/3]; % Select point in center
   rc=[];
   for ff=1:Nb  % Find the position and the normal vector at the center of body's elements
      elknxyzb=nodestmp(topologytmp(ff,:),:);
      [psi, xq, yq, zq, nx, ny, nz]=elemshapestl([elknxyzb ones(nknel*(isquad+1),1)],IP);
      rc=[rc ; xq yq zq nx ny nz]; % Save centers of elements and normal vectors for later use
   end
   
   centr=mean(nodestmp); % mass center of the body
   % find the farthest element
   [dummy,elmax]=max(sum((rc(:,1:3)-ones(Nb,1)*centr).^2,2));
   rcm=rc(elmax(1),:);
   %Cdata(elmax,:)=Cdata(elmax,:); %%%%%%%%%%
   Cdata(elmax,:)=Cdata(elmax,:)*0.6;%[Cdata(elmax,2) Cdata(elmax,1) Cdata(elmax,3)]; %indicate the chosen element by a color change
   ferel=topologytmp(elmax(1),:);
   
   % Are the radiusvector and the farthest element's normal vector in the same hemi-space?
   direction=sign((rcm(1)-centr(1))*rcm(4)+...
                  (rcm(2)-centr(2))*rcm(5)+...
                  (rcm(3)-centr(3))*rcm(6));
               
   if direction < 0
      % change all normal vectors in the body if they point inwards               
      topology(tx,1:nknel)=fliplr(topology(tx,1:nknel));
      topologytmp(:,1:nknel)=fliplr(topologytmp(:,1:nknel));
      %Change also the list of open segments to match the vector flip
      segmcorr1=find(segmopen(:,2)==nknel-1 & segmopen(:,3)==bb);
      segmcorr2=find(segmopen(:,2)==1 & segmopen(:,3)==bb);
      segmopen(segmcorr1,2)=1;
      segmopen(segmcorr2,2)=nknel-1;
      %Change also the list of segments with more than 2 elements attached
      %to match the vector flip (VCH 08/2008)
      segmcorr1=find(segmextraele(:,2)==nknel-1 & segmextraele(:,3)==bb);
      segmcorr2=find(segmextraele(:,2)==1 & segmextraele(:,3)==bb);
      segmextraele(segmcorr1,2)=1;
      segmextraele(segmcorr2,2)=nknel-1;
   end
   
   figure; % one figure for each body, in case they hide one another
   patch('faces',topologytmp(:,1:nknel),'vertices',nodestmp,'FaceVertexCData',Cdata,'FaceColor','flat');
   view(-135,45);
%   axis equal;
   rotate3d on;
   hold on;
   title(['Body nr. ' num2str(bb)]); %%%%%%%%%%
   xlabel('x');ylabel('y');zlabel('z');
   
   hold on;
   if isquad==1
      temp=topologytmp(:,5:8);
      plot3(nodestmp(temp(1:end),1),nodestmp(temp(1:end),2),nodestmp(temp(1:end),3),'k.');
   end
   
   % draw normal vector on the farthest element (debugging)
%   h=line([rcm(1);rcm(4)+rcm(1)],[rcm(2);rcm(5)+rcm(2)],[rcm(3);rcm(6)+rcm(3)]);
%   set(h,'LineWidth',1,'Color',[0 0 1]);
   if ~isempty(segmopen) % draw lines where the body is open
      openbod=segmopen(find(segmopen(:,end)==bb),1:2);
      for ee=1:length(openbod)
         segtmp=topologytmp(find(tx==openbod(ee,1)),mod((openbod(ee,2)-1):(openbod(ee,2)),nknel)+1);
         h=line(nodestmp(segtmp,1),nodestmp(segtmp,2),nodestmp(segtmp,3));
%         set(h,'LineWidth',2,'Color',[0 0 0]); %%%%%%%%%%
         set(h,'LineWidth',1.5,'Color',[1 0 0]);
      end
   end

   if ~isempty(segmextraele) % draw lines where more than two elements are attached
      openbod=segmextraele(find(segmextraele(:,end)==bb),1:2);
      for ee=1:length(openbod)
         segtmp=topologytmp(find(tx==openbod(ee,1)),mod((openbod(ee,2)-1):(openbod(ee,2)),nknel)+1);
         h=line(nodestmp(segtmp,1),nodestmp(segtmp,2),nodestmp(segtmp,3));
         set(h,'LineWidth',2,'Color',[0 1 1]);
      end
   end

end

topologyb=topology;
nodesb=nodes;

toposhrinkb=[];
for rr=1:N
   for cc=1:nknel
      col=mod((cc-1):cc,nknel)+1;
      if topologyb(rr,col(1))==topologyb(rr,col(2))
         topotmp=topologyb(rr,:);
         topotmp(cc)=[];
         toposhrinkb=[toposhrinkb ; topotmp];
      end
   end
end

tim=toc;



function [segmout,topology]=elstep(topology,segmin,body,nknel,N,isquad)
% find element number and position of the segment that matches a given segment
%    -topology: numbers of the nodes in each element, ordered.
%               One element per row. This function re-orders them if necessary
%    -segmin:   row vector with input segment description:
%               element number and position within the element.
%    -body:     current body number.
%    -nknel:    number of nodes per element.
%    -N:        number of elements.
%
%    -segmout:  row vector with output segment description:
%               element number and position within the element.
%               It is [0 1] if an open edge is found, [1 0] if the
%               input segment is duplicated and negative if the
%               new element direction has been corrected.

segbad=[];
segood=[];
segmout=[];

% find all attached segments, matching or not the original segment's direction
for kk=1:nknel
    if isquad==1
        col=mod((kk-1):(nknel+kk-2),nknel)+1;

        tmp1=ones(N,1)*topology(segmin(1),[mod((segmin(2)-1),nknel) mod((segmin(2)-1),nknel)+4 mod(segmin(2),nknel)]+1);
        tmp1(segmin(1),:)=[0 0 0];
        tmp_g=find(sum((abs(topology(:,[col(2) col(1)+4 col(1)])-tmp1))')==0);
        tmp_b=find(sum((abs(topology(:,[col(1) col(1)+4 col(2)])-tmp1))')==0);

        segood=[segood ; tmp_g'  kk*ones(length(tmp_g),1)];
        segbad=[segbad ; tmp_b'  kk*ones(length(tmp_b),1)];
    else
        col=mod((kk-1):(nknel+kk-2),nknel)+1;

        tmp1=ones(N,1)*topology(segmin(1),mod((segmin(2)-1):segmin(2),nknel)+1);
        tmp1(segmin(1),:)=[0 0];
        tmp_g=find(sum((abs(topology(:,col(2:-1:1))-tmp1))')==0);
        tmp_b=find(sum((abs(topology(:,col(1:2))-tmp1))')==0);

        segood=[segood ; tmp_g'  kk*ones(length(tmp_g),1)];
        segbad=[segbad ; tmp_b'  kk*ones(length(tmp_b),1)];
    end
end

% consider all options and correct direction if necessary
szgood=size(segood,1);
szbad=size(segbad,1);
if (szgood+szbad)>1
    if topology(segmin(1),mod((segmin(2)-1),nknel)+1)==topology(segmin(1),mod(segmin(2),nknel)+1)
        segmout=[1 0];
    else  % More than two adjoining elements

        % START: newly added code to cope with more than two adjoining elements (VCH 08/2008)
        if szgood>=1
            segmout=segood;
        end
        if szbad>=1
            %note: this specializes to triangular elements
%             for ii=1:size(segbad,1)
%                 dummy=topology(segbad(ii,1),:);
%                 col=mod((segbad(ii,2)-1):(nknel+segbad(ii,2)-2),nknel)+1;
%                 topology(segbad(ii,1),col(1))=dummy(col(2));
%                 topology(segbad(ii,1),col(2))=dummy(col(1));
%                 topology(segbad(ii,1),col(3))=dummy(col(3));
%                 %topology(segbad(ii,1),col(4))=dummy(col(3));
%                 if isquad==1
%                     topology(segbad(ii,1),col(2)+4)=dummy(col(4)+4);
%                     topology(segbad(ii,1),col(4)+4)=dummy(col(2)+4);
%                 end
%                 if segbad(ii,2)==nknel-1, segbad(ii,2)=1; end     % VCH 08/2008  ?????????????????
%                 if segbad(ii,2)==1, segbad(ii,2)=nknel-1; end     % VCH 08/2008
%                 segmout=[segmout ;-segbad(ii,:)];
%             end
            segmout=[segmout ;segbad];
        end
        % END: newly added code to cope with more than two adjoining elements (VCH 08/2008)

    end
elseif (szgood+szbad)==0
    if topology(segmin(1),mod((segmin(2)-1),nknel)+1)==topology(segmin(1),mod(segmin(2),nknel)+1)
        segmout=[1 0];
    else
        segmout=[0 1];
    end
elseif szgood==1
    segmout=segood;
elseif szbad==1
%     %note: this specializes to triangular elements
%     dummy=topology(segbad(1),:);
%     col=mod((segbad(2)-1):(nknel+segbad(2)-2),nknel)+1;
%     topology(segbad(1),col(1))=dummy(col(2));
%     topology(segbad(1),col(2))=dummy(col(1));
%     topology(segbad(1),col(3))=dummy(col(3));
%     %topology(segbad(1),col(4))=dummy(col(3));
%     if isquad==1
%         topology(segbad(1),col(2)+4)=dummy(col(4)+4);
%         topology(segbad(1),col(4)+4)=dummy(col(2)+4);
%     end
%     if segbad(2)==nknel-1, segbad(2)=1; end     % VCH 08/2008 ????????????????????????
%     if segbad(2)==1, segbad(2)=nknel-1; end     % VCH 08/2008
%     segmout=-segbad;
    segmout=segbad;
end

function [psi, xq, yq, zq, nx, ny, nz]=elemshapestl(elknxyzb,IP)
% Calculate coordinates and normalvector at 'x'
% in an isoparametric triangular element. (linear).
%
%  Input:
%    IP:     real matrix containing the ordinates for the
%             integration points (IP's) x values in 1st column y values in 2nd
%    elknxyzb: real matrix, one row for each node in the element
%             each row contains (x,y,z,body) for the node
%
%  Output:
%    psi:     Real matrix containing one shape function in each
%             column (no. of columns=no. of nodes in an element)
%    xq:      Real column vector containing the x-coordinates
%             of the element for each IP
%    yq:      Real column vector containing the y-coordinates
%             of the element for each IP
%    zq:      Real column vector containing the z-coordinates of
%             the element for each IP
%    nx,ny,nz: 
%             The column vectors containing respectively the x,y
%             and z-coordinates of the normal vector.
%             Note that the normal vector does not have unit length.
%             This is done to avoid division by the jacobian which
%             occasionally is zero. Subsequent multiplication must
%             therefore be omitted.

[nknel,dummy]=size(elknxyzb);
if nknel==3
   % Linear Shape functions:
   psi=[IP(:,1) ...
        IP(:,2) ...
        1-IP(:,1)-IP(:,2)];
   % Linear Shape function derivatives
   % Well, they are actually constant
   ettaller=ones(size(IP,1),1);
   nuller=zeros(size(IP,1),1);
   
   dNx=[ettaller nuller -ettaller];
   dNy=[nuller ettaller -ettaller];
   
%    dNx(:,1)=ones(size(IP,1),1);
%    dNx(:,2)=zeros(size(IP,1),1);
%    dNx(:,3)=-ones(size(IP,1),1);
%    
%    dNy(:,1)=zeros(size(IP,1),1);
%    dNy(:,2)=ones(size(IP,1),1);
%    dNy(:,3)=-ones(size(IP,1),1);
  else
   error('Only linear shape functions are implemented');
end

xq=psi*elknxyzb(:,1);
yq=psi*elknxyzb(:,2);
zq=psi*elknxyzb(:,3);

% Elements of the Jacobian matrix for all input points: dr/de1 and dr/de2, r=(x,y,z).
dxde1=dNx*elknxyzb(:,1);
dyde1=dNx*elknxyzb(:,2);
dzde1=dNx*elknxyzb(:,3);

dxde2=dNy*elknxyzb(:,1);
dyde2=dNy*elknxyzb(:,2);
dzde2=dNy*elknxyzb(:,3);

% Normal vector at input points: scalar product (dr/de1 x dr/de2).
IntExt=sign(elknxyzb(1,4)); % Change sign if the domain is interior. (VC)
nx=(dyde1.*dzde2-dzde1.*dyde2)*IntExt;
ny=(dzde1.*dxde2-dxde1.*dzde2)*IntExt;
nz=(dxde1.*dyde2-dyde1.*dxde2)*IntExt;


