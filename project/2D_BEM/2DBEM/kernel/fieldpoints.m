function [Ap,Bp,CConst]=fieldpoints(xyb,topology,k,betaP,xy);

% [Ap,Bp,CConst]=fieldpoints(xyb,topology,k,sigma,xy);
%
% Calculates rows of coefficients for a set of field points.
% It admits a plane with finite or infinite impedance.
%
% Input variables:
%   -xyb:     node positions, first column is the x-coordinate, second
%             column is y-coordinate, and third column is the body
%             number to which the node belongs to.
%   -topology:each row contains the node numbers (row number in xyb) of the
%             nodes in one element. The last column is the body number the
%             element belongs to.
%   -k:       wavenumber.
%   -betaP:   normalised admittance of the plane, at k. If its value is 0,
%             a rigid plane is considered (infinite impedance); if its value
%             is NaN (not-a-number), free-field is assumed.
%   -xy:      field points to calculate. Column 1, x coordinate; column 2,
%             y coordinate.
%
% Output variable:
%   -Ap:      coefficient matrix for pressure.
%   -Bp:      coefficient matrix for normal velocity.
%   -CConst:  C constants, in this case they are not substracted from the matrix. 
%
% It is possible to calculate for the interior or exterior domain, or a
% combination ob both. The sign of the body numbers indicates whether the
% interior or exterior domain to that body must be considered.
%
% When non-continuous boundary conditions exist, the columns in
% matrix B can be expanded. An imaginary part must be added to
% the body numbers in 'topology' to indicate what boundary area
% the element belongs to (1,2,3...). The function 'bound2D' can be
% used to set up the equations and define the boundary areas.(VC)

% Vicente Cutanda Henriquez 5-2001.
% VCH 4-2017: B-matrix expansion, BC splitting. 

[M,ncolxyb] = size(xyb);
[N,nknel] = size(topology);
np=size(xy,1);
nknel=nknel-1;
xy=[xy ones(np,1)];

NumBodies = max(abs(xyb(:,3)));
if NumBodies ~= max(abs(real(topology(:,nknel+1))))
   error('Error: The input arrays must include consistent body numbers - Calculation aborted');
end

% Check geometry files to find dummy elements and nodes
if isnan(betaP)
   operate='n';
else
   operate='y';
end
[xybnd,topologynd,xydum,topodum,xynodum]=nodummy(xyb,topology,operate);
lndu=size(xybnd,1);
ledu=size(topologynd,1);

% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
CConst=2*pi*ones(np,1);
Ap=zeros(np,lndu);
if isreal(topology) | operate=='y' % calculations with impedance plane are not combined with B matrix extension
    BCtopo=topologynd;
    Bp=Ap;
else 
    BCelem=imag(topology(:,end));
    topologynd=real(topology);
    topology=real(topology);
    % see help in function 'bound' (VC)
    [BCtopo,BCnodeA,BCnodeB]=bound2D(xyb,topology,BCelem,'n');
    Bp=zeros(np,size(BCnodeB,1));
end

% disp([' 2D BEM field point calculation, ' num2str(np) ' points, k = ' num2str(k)])

% coefficient matrix calculation
for jj=1:np
   for nde=1:N
      elknxyb=xyb(topology(nde,1:nknel),:);
      [F1A,F1B,CK]=intF1(xy(jj,:),elknxyb,k,betaP);
      CConst(jj)=CConst(jj)+CK;
      if topodum(nde)~=0
         Ap(jj,topologynd(topodum(nde),1:nknel))=Ap(jj,topologynd(topodum(nde),1:nknel))+F1A(1:nknel);
         Bp(jj,BCtopo(topodum(nde),1:nknel))=Bp(jj,BCtopo(topodum(nde),1:nknel))+F1B(1:nknel);
      end
   end
%    if np>1 && jj/100==round(jj/100)
%        disp([' 2D BEM field point calculation, point ' num2str(jj) ' of ' num2str(np) ...
%          ', C constant = ' num2str(CConst(jj)/pi) ' * pi'])
%    end
end
