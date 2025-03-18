function [bp,wf]=nsingrule(n,elknrzb,przb)

% function [bp,wf]=nsingrule(n,elknrzb,przb)
%
% 1-D integration rule for near singular integrals.
% Returns the weights and base points for integration with
% a near singularity. Recursive version.
% Recalculates the distance of the subelements to every subelement
%
% Input:
%   -n:       Order of Gauss-Legendre rule in each subdivision.
%   -elknrzb: real matrix, one row for each node in the element
%             each row contains (rho,z,body) for the node
%   -przb:    close point coordinates (rho,z,body)
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

[nknel,dummy]=size(elknrzb);

% Uses Gauss-Legendre in each direction for each subdivision
[bpT,wfT]=gaussrule(n);

% Call to recursive subdivision function
DivsIN=[-1 1]; % limits of the whole element
if nknel==3 % Quadratic elements
   % Calculate length of the subelement. Approximate expression for the
   % length of a circular arc: l=(8*b-a)/3, with a:chord, b:chord of half
   % the arc. (From Bronshtein-Semendiaev Mathematics manual, p 140)
    a=abs((elknrzb(1,1)+j*elknrzb(1,2))-(elknrzb(3,1)+j*elknrzb(3,2)));
    b=abs((elknrzb(1,1)+j*elknrzb(1,2))-(elknrzb(2,1)+j*elknrzb(2,2)));
    lengthELE=(8*b-a)/3;
elseif nknel==2 % linear element, no arc to calculate
    lengthELE=abs((elknrzb(1,1)+j*elknrzb(1,2))-(elknrzb(2,1)+j*elknrzb(2,2)));
end
DivsOUT=subdivideAXI(DivsIN,elknrzb,przb,lengthELE);

Ndivs=size(DivsOUT,1);
if Ndivs>3 & min(sqrt((elknrzb(:,1)-przb(1)).^2+(elknrzb(:,2)-przb(2)).^2))>eps
%     disp(['Near Singular Integration Performed, nr. divisions: ' num2str(Ndivs)]);
end

% Fill subelements with integration points
bp=[];wf=[];
for ii=1:size(DivsOUT,1)
   Dim=DivsOUT(ii,2)-DivsOUT(ii,1);
   Shift=(DivsOUT(ii,1)+DivsOUT(ii,2))/2;
   bp=[bp ; bpT/2*Dim+Shift];
   wf=[wf ; wfT/2*Dim];
end



function DivsOUT=subdivideAXI(DivsIN,elknrzb,przb,lengthELE)

% DivsOUT=subdivideAXI(DivsIN,elknrzb,przb,lengthELE)
%
% Recursive function performing element subdivision.
%
% Input:
%  -DivsIN:   Contains the coordinates of left and right limits of subelements
%             arranged with one subelement per row.
%   -elknrzb: real matrix, one row for each node in the element
%             each row contains (rho,z,body) for the node
%   -przb:    close point coordinates (rho,z,body)
%   -lengthELE: approximate length of the whole element
%
% Output:
%  -DivsOUT:  The list of subelements created from DivsIN, arranged as DivsIN.

% Vicente Cutanda Henriquez 11-2010
[nknel,dummy]=size(elknrzb);

DivsOUT=[];
for ii=1:size(DivsIN,1) % loop over all current subelements
   % The two vertices of the subelement and the central point
   P1=DivsIN(ii,1); P2=DivsIN(ii,2);
   Pcent=(P1+P2)/2;
   % Convert to global coordinates
   [psi, rhoq, zq]=elemshape(elknrzb,[P1 Pcent P2]);
   % Calculate distance from element centre to calculation point
   DistSUBELE_NS=sqrt((rhoq(2)-przb(1))^2 + (zq(2)-przb(2))^2);
if nknel==3 % Quadratic elements
   % Calculate length of the subelement. Approximate expression for the
   % length of a circular arc: l=(8*b-a)/3, with a:chord, b:chord of half
   % the arc. (From Bronshtein-Semendiaev Mathematics manual, p 140)
   a=abs((rhoq(1)+j*zq(1))-(rhoq(3)+j*zq(3)));
   b=abs((rhoq(1)+j*zq(1))-(rhoq(2)+j*zq(2)));
   lengthSUBELE=(8*b-a)/3;
elseif nknel==2 % linear element, no arc to calculate
   lengthSUBELE=abs((rhoq(1)+j*zq(1))-(rhoq(2)+j*zq(2)));
end
   % Limit the subdivision if the calculation point is too close 
   % or is a collocation point (singular integration):
   if lengthSUBELE<lengthELE*1e-3, lengthSUBELE=0; end
   NewDivs=[];
   if lengthSUBELE>DistSUBELE_NS/2 % Condition for recursion
      % Subdivide in two a particular subelement
      NewDivs=[P1 Pcent;Pcent P2];
      % Recursion call
      NewDivs=subdivideAXI(NewDivs,elknrzb,przb,lengthELE);
      DivsOUT=[DivsOUT ; NewDivs];
   else
      DivsOUT=[DivsOUT ; DivsIN(ii,:)];
   end
end


