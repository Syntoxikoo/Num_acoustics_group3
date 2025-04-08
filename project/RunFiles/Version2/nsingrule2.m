function [bp,wf]=nsingrule2(n,xising,ndist)

% function [bp,wf]=nsingrule2(n,xising,ndist)
%
% 1-D integration rule for near singular integrals.
% Returns the weights and base points for integration with
% a near singularity. Recursive version.
%
% Input:
%   -n:        Order of Gauss-Legendre rule in each subdivision.
%   -xising:   Local coordinate of the projection of the near point.
%   -ndist:    Normalized (to element size) distance of near point to the element.
%
% Output:
%   -bp, wf:   Base points and weights for the integration.
%
% Details on this technique are given in:
%
% "Numerical transducer modelling", by Vicente Cutanda
% Proceedings of the 6th International Congress on Sound and Vibration
% Copenhaguen, July 1999
%
% "On the modeling of narrow gaps using the standard boundary element method"
% Vicente Cutanda, Peter M. Juhl and Finn Jacobsen
% J. Acoust. Soc. Am. 109 (4), April 2001, 1296-1303.
%
% “Acoustic boundary element method formulation with treatment of nearly
% singular integrands by element subdivision”
% Vicente Cutanda Henríquez and Peter Juhl
% 19th International Congress on Acoustics, Madrid, 2-7 September 2007

% Version by Vicente Cutanda Henriquez 11-2010

% This function is not used at the moment. It would need a different implementation
% where the distance and projection of the near-singular point are calculated 


MinSize=ndist;    % Depth of subdivision (see subfunction)

MaxDist=1;      % Subdivision distance threshold (see subfunction)

% Uses Gauss-Legendre in each direction for each subdivision
[bpT,wfT]=gaussrule(n);

% Call to recursive subdivision function
DivsIN=[-1 1]; % limits of the whole element
DivsOUT=subdivideAXI(DivsIN,xising,MinSize,MaxDist);

% Fill subelements with integration points
bp=[];wf=[];
for ii=1:size(DivsOUT,1)
   Dim=DivsOUT(ii,2)-DivsOUT(ii,1);
   Shift=(DivsOUT(ii,1)+DivsOUT(ii,2))/2;
   bp=[bp ; bpT/2*Dim+Shift];
   wf=[wf ; wfT/2*Dim];
end



function DivsOUT=subdivideAXI(DivsIN,xising,MinSize,MaxDist)

% DivsOUT=subdivideAXI(DivsIN,xising,MinSize,MaxDist)
%
% Recursive function performing element subdivision.
%
% Input:
%  -DivsIN:  Contains the coordinates of left and right limits of subelements arranged
%            with one subelement per row.
%  -xising:  Local coordinates of the projection of the near point.
%  -MinSize: Minimum legth of the subelements, to put an end to the recursion.
%  -MaxDist: If a subelement's center is less than MaxDist times its side away
%            from the projection of the near point and its area if larger than
%            MinSize, it is subdivided.
%
% Output:
%  -DivsOUT: The list of subelements created from DivsIN, arranged as DivsIN.

% Vicente Cutanda Henriquez 11-2010

DivsOUT=[];
for ii=1:size(DivsIN,1) % loop over all current subelements
   % The two vertices of the subelement
   P1=DivsIN(ii,1); P2=DivsIN(ii,2);
   % Calculate subelement dimensions
   Pcent=(P1+P2)/2; lengthSUBELE=abs(P2-P1);
   % Calculate distance from element centre to near-singular point projection
   DistSUBELE_NS=min(abs(DivsIN(ii,:)-xising));
   NewDivs=[];
   if lengthSUBELE>MinSize & DistSUBELE_NS<MaxDist*lengthSUBELE % Conditions for recursion
      % Subdivide in two a particular subelement
      NewDivs=[P1 Pcent;Pcent P2];
      % Recursion call
      NewDivs=subdivideAXI(NewDivs,xising,MinSize,MaxDist);
      DivsOUT=[DivsOUT ; NewDivs];
   else
      DivsOUT=[DivsOUT ; DivsIN(ii,:)];
   end
end
