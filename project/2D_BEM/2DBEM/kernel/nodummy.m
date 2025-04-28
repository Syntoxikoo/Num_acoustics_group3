function [xybnd,topologynd,xydum,topodum,xynodum]=nodummy(xyb,topology,operate);

% [xybnd,topologynd,xydum,topodum,xynodum]=nodummy(xyb,topology,operate);
% 
% Reduces the geometry matrices 'xyb' and 'topology' by eliminating dummy elements (all y = 0).
% Input:
%     -operate:if 'y', the function works, otherwise it does nothing but fill the output variables
%              with default values.
% Output:
%     -xydum:  vector with the same rows as 'xybnd', its elements are 1 if the corresponding node 
%              is on the y=0 plane, 0 otherwise.
%     -topodum:vector with the same rows as 'topology'. If the corresponding element is not dummy,
%              it gives its row number in 'topologynd'.
%     -xynodum:indexes in xyb which have a corresponding no dummy node.

% Vicente Cutanda 5-2001.


if operate=='y'
   xybnd=[];
   topologynd=[];
   topodum=zeros(size(topology,1),1);
   xydum=[];
   xynodum=[];
   tk=zeros(size(xyb,1),1);
   for nn=1:size(topology,1)
       if any(xyb(topology(nn,1:end-1),2)>eps) %| xyb(topology(nn,1),1)<xyb(topology(nn,2),1)
         topologynd=[topologynd ; topology(nn,:)];
         topodum(nn)=size(topologynd,1);
         for bb=1:length(topology(nn,1:end-1))
            if tk(topology(nn,bb))==0
               xybnd=[xybnd ; xyb(topology(nn,bb),:)];
               xynodum=[xynodum ; topology(nn,bb)];
               if xyb(topology(nn,bb),2)<eps, aa=1; else, aa=0; end
               xydum=[xydum; aa];
               tk(topology(nn,bb))=size(xybnd,1);
            else
               if xyb(topology(nn,bb),2)<eps
                  xydum(tk(topology(nn,bb)))=1;
               end
            end
            topologynd(end,bb)=tk(topology(nn,bb));
         end
      end
   end
else
   xybnd=xyb;
   xynodum=(1:size(xyb,1))';
   topologynd=topology;
   xydum=zeros(size(xyb,1),1);
   topodum=[1:size(topology,1)]';
end
