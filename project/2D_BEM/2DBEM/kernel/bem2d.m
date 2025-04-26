function [A,B]=bem2d(xyb,topology,k,betaP,varargin);

% [A,B]=bem2d(xyb,topology,k,betaP{,chiefpoints});
%
% Calculates the coefficient matrix for the 2D BEM formulation.
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
%   -chiefpoints: Like 'xyb', but contains CHIEF points instead
%             One row for each chief point
%
% Output variable:
%   -A:       coefficient matrix for the pressure.
%   -B:       coefficient matrix for the normal velocity.
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
% Vicente Cutanda Henriquez 03-2011, version including CHIEF points
% VCH 4-2017: B-matrix expansion, BC splitting. 


[M,ncolxyb] = size(xyb);
[N,nknel] = size(topology);
nknel=nknel-1;
if ncolxyb<3
   error('Error: The input arrays must include body numbers - Calculation aborted');
end

if nargin<5
   chiefpoints=[];
else
   chiefpoints=varargin{1};
end
nchiefp=size(chiefpoints,1);

NumBodies = max(abs(xyb(:,3)));
if NumBodies ~= max(abs(real(topology(:,nknel+1))));
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
% ledu=size(topologynd,1);


% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
CConst=2*pi*(1+sign(xybnd(:,3)))/2; % 2*pi or 0 for exterior/interior domain.
A=zeros(lndu+nchiefp,lndu);
if isreal(topology) | operate=='y' % calculations with impedance plane are not combined with B matrix extension
    BCtopo=topologynd;
    B=A;
else 
    BCelem=imag(topology(:,end));
    topologynd=real(topology);
    topology=real(topology);
    % see help in function 'bound' (VC)
    [BCtopo,BCnodeA,BCnodeB]=bound2D(xyb,topology,BCelem,'n');
    B=zeros(lndu+nchiefp,size(BCnodeB,1));
end

disp([' 2D BEM calculation, ' num2str(lndu+nchiefp) ' points, k = ' num2str(k)])

% coefficient matrix calculation
for jj=1:lndu+nchiefp
    if jj <= lndu
        pxyb=xybnd(jj,:);
    else
        pxyb=chiefpoints(jj-lndu,:);
    end
    for nde=1:N
        elknxyb=xyb(topology(nde,1:nknel),:);
        [F1A,F1B,CK]=intF1(pxyb,elknxyb,k,betaP);
        if jj <= lndu
            CConst(jj)=CConst(jj)-CK*(1+xydum(jj));% Nodes in contact with the plane have different C.
        end
        if topodum(nde)~=0
            A(jj,topologynd(topodum(nde),1:nknel))=A(jj,topologynd(topodum(nde),1:nknel))+F1A(1:nknel);
            B(jj,BCtopo(topodum(nde),1:nknel))=B(jj,BCtopo(topodum(nde),1:nknel))+F1B(1:nknel);
        end
    end
    if jj <= lndu
%         if jj/100==round(jj/100)
%             disp([' 2D BEM calculation, row ' num2str(jj) ' of ' num2str(lndu+nchiefp) ...
%                 ', C constant = ' num2str(CConst(jj)/pi) ' * pi'])
%         end
        A(jj,jj)=A(jj,jj)-CConst(jj);
    end
end
