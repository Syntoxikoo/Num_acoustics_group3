function is_close=nsingcheck(pxyzb,elknxyzb)

% is_close=nsingcheck(pxyzb,elknxyzb)
%
%  Checks if the calculation point, given in 'pxyzb',
%  is close to the current element, given in 'elknxyzb'.
%  If so, the output variable is set to "1", othervise "0".
%  This very simnplified version is intended as a pre-check to decide if
%  an itegration should be handled as near-singular. It works on triangular
%  and quadrilateral elements

% Vicente Cutanda 12-2010

elknxyz=elknxyzb(:,1:3);pxyz(1,1:3)=pxyzb(1:3);

% Get size of the maximum distance between nodes
maxnoddist=max(sqrt(sum(diff([elknxyz; elknxyz(1:2:end,:); elknxyz(2:2:end,:)]).^2,2)));
% distance from the calculation point to the farthest node in the element:
maxdist=max(sqrt(sum((elknxyz-ones(size(elknxyz,1),1)*pxyz).^2,2)));

%testfactor=1.1;
testfactor=2;
%testfactor=4;

if maxdist>maxnoddist*testfactor % it should not get nodes from adjoining elements
   is_close=0;
else
   is_close=1;
end

