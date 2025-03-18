function [topology]=rearrange(topology,n);
% [topology]=rearrange(topology,n);
% rearrange the nodal order in element(s) n to shift normal direction

[nel,nknel]=size(topology);
if nknel~=4 & nknel~=8
   error('input should be linear or quadric quadrilaterals');
end

for jj=1:length(n)
   oldelement(1,:)=topology(n(jj),:);
   topology(n(jj),2)=oldelement(1,4);
   topology(n(jj),4)=oldelement(1,2);
   if nknel==8
      topology(n(jj),5)=oldelement(1,8);
      topology(n(jj),6)=oldelement(1,7);
      topology(n(jj),7)=oldelement(1,6);
      topology(n(jj),8)=oldelement(1,5);
   end
   
   
end
