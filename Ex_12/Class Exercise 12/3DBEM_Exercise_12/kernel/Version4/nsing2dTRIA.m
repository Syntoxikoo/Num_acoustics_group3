function nsingIP=nsing2dTRIA(n,pxyzb,elknxyzb,varargin)

% nsingIP=nsing2dTRIA(n,pxyzb,elknxyzb,Tole);
%
% 2-D integration rule for near singular integrals. Version for triangular elements.
%
% Input:
%   -n:        Order of Gauss-Legendre rule in each subdivision (3, 4, 7).
%    pxyzb:    real vector containing the (x,y,z,body) values for
%              the point 'P'.
%    elknxyzb: real matrix, one row for each node in the element
%              each row contains (x,y,z,body) for the node.
%    -Tole:    Minimum relative size of the subelements for near-singular
%              integrals. Default:1e-6. If Tole=0, no subdivision is made.
%
% Output:
%   -nsingIP:  Integration points, 1st column e1-values; 2nd column e2-values; 3rd column weigths.

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

% Version by Vicente Cutanda Henriquez 12-2010

% VCH 03-2015: filling loop modified for higher speed.

if nargin==3
   Tole=1e-6;
else
   Tole=varargin{1};
end

see=0;            % Produce a 3D figure of the integration points (debugging)

% Generic subelement mesh:
IPgen=gaussTRI(n);

if Tole==0 % exit with undivided element
    nsingIP=IPgen;
    return
end

% Get the maximum distance between nodes as a reference:
sizeELE=max(sqrt(sum(diff([elknxyzb(:,1:3); elknxyzb(1,1:3)]).^2,2)));

% Call to recursive subdivision function
DivsIN=[1 0 0 1 0 0]; % The three nodes of the full element (local coordinates)
DivsOUT=subdivideTRIA(DivsIN,pxyzb,elknxyzb,sizeELE,Tole);

% Fill subelements with integration points
nsingIP=zeros(size(DivsOUT,1)*size(IPgen,1),3);
for ii=1:size(DivsOUT,1) % Scale and displace the IPs and weights for every subelement
   Dim=min(sqrt(sum(diff([DivsOUT(ii,1:2);DivsOUT(ii,3:4);DivsOUT(ii,5:6);DivsOUT(ii,1:2)]).^2,2)));
   Pcent=(DivsOUT(ii,1:2)+DivsOUT(ii,3:4)+DivsOUT(ii,5:6))/3;
   nsingIP(((ii-1)*size(IPgen,1)+1):ii*size(IPgen,1),:)=...
           [(IPgen(:,1)-1/3)*Dim+Pcent(1) (IPgen(:,2)-1/3)*Dim+Pcent(2) IPgen(:,3)*Dim^2];
end

% if size(DivsOUT,1)>1
%     disp(['Near-singular integration - Subelements: ' num2str(size(DivsOUT,1)) '    IPs: ' num2str(size(nsingIP,1))]);
% end

% To see the result (debugging):
if see
   figure;
   plot3(nsingIP(:,1),nsingIP(:,2),nsingIP(:,3)*30,'b.',[1;0;0],[0;1;0],[0;0;0],'r*',...
      [DivsOUT(:,1);DivsOUT(:,3);DivsOUT(:,5)],[DivsOUT(:,2);DivsOUT(:,4);DivsOUT(:,6)],zeros(3*size(DivsOUT,1),1),'g*');
   grid; view(0,90);
   title(['Subelements: ' num2str(size(DivsOUT,1)) '    IPs: ' num2str(size(nsingIP,1))])
   rotate3d on
end


function DivsOUT=subdivideTRIA(DivsIN,pxyzb,elknxyzb,sizeELE,Tole)

% DivsOUT=subdivideTRIA(DivsIN,pxyzb,elknxyzb,sizeELE,Tole)
%
% Recursive function performing element subdivision. Works on triangular elements.
%
% Input:
%  -DivsIN:  Contains the (2) coordinates the 3 vertices of subelements arranged
%            with one subelement per row:
%            [x1 y1  x2 y2 ; ...]
%            where the vertices are numbered clockwise.
%    pxyzb:    real vector containing the (x,y,z,body) values for
%              the point 'P'.
%    elknxyzb: real matrix, one row for each node in the element
%              each row contains (x,y,z,body) for the node.
%    -Tole:    Minimum relative size of the subelements for near-singular
%              integrals. Default:1e-6.
%
% Output:
%  -DivsOUT: The list of subelements created from DivsIN, arranged as DivsIN.

% Vicente Cutanda 12-2010

DivsOUT=[];
for ii=1:size(DivsIN,1) % loop over all current subelements
   % Two of the three vertices of the subelement
   P1xy=DivsIN(ii,1:2); P2xy=DivsIN(ii,3:4);  P3xy=DivsIN(ii,5:6);
   % New points to get dimensions and subdivide
   P12xy=(P1xy+P2xy)/2; P23xy=(P2xy+P3xy)/2; P31xy=(P3xy+P1xy)/2;
   Pcentxy=(P1xy+P2xy+P3xy)/3;
   [psi, xq, yq, zq]=elemshapetri(elknxyzb,[Pcentxy; P1xy; P2xy; P3xy]);
   % Calculate minimum distance from element nodes and centre to calculation point
   DistSUBELE_NS=min(sqrt((xq-pxyzb(1)).^2 + (yq-pxyzb(2)).^2 + (zq-pxyzb(3)).^2));
   % Calculate size of the subelement, as the maximum distance between nodes
   subelknxyz=[xq(2:end,1) yq(2:end,1) zq(2:end,1)];
   sizeSUBELE=max(sqrt(sum(diff([subelknxyz; subelknxyz(1,:)]).^2,2)));

   % Limit the subdivision if the calculation point is too close 
   if sizeSUBELE<sizeELE*Tole, sizeSUBELE=0; end
   
   NewDivs=[];
   Ftest=4;   %  Test constant to increase the division on subelements away from the NS point (initially 1) VCH 3-10-2012
   if sizeSUBELE>DistSUBELE_NS/Ftest % Condition for recursion <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      % Subdivide in four a particular subelement
      NewDivs=[P1xy P12xy P31xy ; P12xy P2xy P23xy ; P31xy P23xy P3xy ; P12xy P23xy P31xy];
      % Recursion call
      NewDivs=subdivideTRIA(NewDivs,pxyzb,elknxyzb,sizeELE,Tole);
      DivsOUT=[DivsOUT ; NewDivs];
   else
      DivsOUT=[DivsOUT ; DivsIN(ii,:)];
   end
end
