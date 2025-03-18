function [xyzbnd,topologynd,xyzdum,topodum,xyznodum]=nodummy3D(xyzb,topology,operate);

% [xyzbnd,topologynd,xyzdum,topodum,xyznodum]=nodummy3D(xyzb,topology,operate);
% 
% Reduces the geometry matrices 'xyzb' and 'topology' by eliminating dummy elements (all z = 0).
% Input:
%     -operate:if 1, the function works, otherwise it does nothing but fill the output variables
%              with default values.
% Output:
%     -xyzdum: vector with the same rows as 'xyzbnd', its elements are 1 if the corresponding node 
%              is on the z=0 plane, 0 otherwise.
%     -topodum:vector with the same rows as 'topology'. If the corresponding element is not dummy,
%              it gives its row number in 'topologynd'.
%     -xyznodum:indexes in xyzb which have a corresponding no dummy node.

% Vicente Cutanda Henriquez 6-2016.


if operate
   xyzbnd=[];
   topologynd=[];
   topodum=zeros(size(topology,1),1);
   xyzdum=[];
%   xyznodum=[];
   xyznodum=zeros(size(xyzb,1),1);
   tk=zeros(size(xyzb,1),1);
   for nn=1:size(topology,1)
       if any(xyzb(topology(nn,1:end-1),3)>eps) %| xyb(topology(nn,1),1)<xyb(topology(nn,2),1)
         topologynd=[topologynd ; topology(nn,:)];
         topodum(nn)=size(topologynd,1);
         for bb=1:length(topology(nn,1:end-1))
            if tk(topology(nn,bb))==0
               xyzbnd=[xyzbnd ; xyzb(topology(nn,bb),:)];
               xyznodum(topology(nn,bb),1)=size(xyzbnd,1);
%                xyznodum=[xyznodum ; topology(nn,bb)];
               if xyzb(topology(nn,bb),3)<eps
                   aa=1;
               else
                   aa=0;
               end
               xyzdum=[xyzdum; aa];
               tk(topology(nn,bb))=size(xyzbnd,1);
            else
               if xyzb(topology(nn,bb),3)<eps
                  xyzdum(tk(topology(nn,bb)))=1;
               end
            end
            topologynd(end,bb)=tk(topology(nn,bb));
         end
      end
   end
else
   xyzbnd=xyzb;
   xyznodum=(1:size(xyzb,1))';
   topologynd=topology;
   xyzdum=zeros(size(xyzb,1),1);
   topodum=[1:size(topology,1)]';
end
