function [A,B]=bem2d(xyb,topology,k);

% [A,B]=bem2d(xyb,topology,k);
%
% Calculates the coefficient matrix for the 2D BEM formulation.
%
% Input variables:
%   -xyb:     node positions, first column is the x-coordinate, second
%             column is y-coordinate, and third column is the body
%             number to which the node belongs to.
%   -topology:each row contains the node numbers (row number in xyb) of the
%             nodes in one element. The last column is the body number the
%             element belongs to.
%   -k:       wavenumber.
%
% Output variable:
%   -A:       coefficient matrix for the pressure.
%   -B:       coefficient matrix for the normal velocity.
%
% It is possible to calculate for the interior or exterior domain, or a
% combination ob both. The sign of the body numbers indicates whether the
% interior or exterior domain to that body must be considered.

% Vicente Cutanda Henriquez 05-2001

[M,ncolxyb] = size(xyb);
[N,nknel] = size(topology);
nknel=nknel-1;
if ncolxyb<3
   error('Error: The input arrays must include body numbers - Calculation aborted');
end

NumBodies = max(abs(xyb(:,3)));
if NumBodies ~= max(abs(topology(:,nknel+1)));
   error('Error: The input arrays must include consistent body numbers - Calculation aborted');
end

% coefficient matrix calculation
CConst=2*pi*(1+sign(xyb(:,3)))/2; % 2*pi or 0 for exterior/interior domain.
%CConst=[CConst ; zeros(nchiefp,1)];
A=zeros(M);
B=zeros(M);
for jj=1:M
    pxyb=xyb(jj,:);
    for nde=1:N
        elknxyb=xyb(topology(nde,1:nknel),:);
        [F1A,F1B,CK]=intF1(pxyb,elknxyb,k);
        CConst(jj)=CConst(jj)-CK;% Nodes in contact with the plane have different C.
        A(jj,topology(nde,1:nknel))=A(jj,topology(nde,1:nknel))+F1A(1:nknel);
        B(jj,topology(nde,1:nknel))=B(jj,topology(nde,1:nknel))+F1B(1:nknel);
    end
    if jj/100==round(jj/100)
        disp([' 2D BEM calculation, row ' num2str(jj) ' of ' num2str(M) ...
            ', C constant = ' num2str(CConst(jj)/pi) ' * pi'])
    end
    A(jj,jj)=A(jj,jj)-CConst(jj);
end
