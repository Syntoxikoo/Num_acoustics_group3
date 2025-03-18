function [Arows,Brows,C]=FieldPntRows(fprz,rzb,Topology,k,m)

%  [Arows,Brows,C]=FieldPnt(fprz,rzb,Topology,k,m)
%
%  Axisymetric BEM calculation.
%  Calculate pressure in field points given a solution on the boundary
%  
%  Input parameters (in SI units):
%  
%    -fprz :    Field points' coordinates. First column are rho-coordinates,
%               second column are z-coordinates.
%    -rzb     : Geometry matrix, one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3 (not necessary if there is just one body)
%    -Topology: Topology matrix
%               One row for each element. Each row contains global node numbers
%               in the column no. corresponding to the local
%               node number. Last column contains the body number.
%    -k       : Wave number (k=2*pi*f/c) (real scalar). 
%    -m       : Circumferential mode number
%
%  Output parameters:
%  
%    -Arows,Brows: Rows of coefficients to apply to the pressure at the 
%                boundary and its normal derivative.
%    -C        : C constants from the calculation. 

% Vicente Cutanda Henríquez, 1-2013.

nfp=size(fprz,1); 
M=size(rzb,1);
[nel, ncols]=size(Topology);
nknel=ncols-1;
cii=0;
Arows=zeros(nfp,M);
Brows=zeros(nfp,M);
C=zeros(nfp,1);
for ff=1:nfp
    A=zeros(1,M);
    B=zeros(1,M);
    cii=4*pi*(1+sign(rzb(1,3)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
    for iel=1:nel
        elknrzb(1:nknel,:)=rzb(Topology(iel,1:nknel),:);
        % Singular part
        [g,h,cjj] = intF2([fprz(ff,:) elknrzb(1,3)],elknrzb,m);
        cii=cii+cjj;
        A(Topology(iel,1:nknel))=A(Topology(iel,1:nknel))+h(1:nknel).';
        B(Topology(iel,1:nknel))=B(Topology(iel,1:nknel))+g(1:nknel).';
        % Oscillating part
        [g,h] = intF1([fprz(ff,:) elknrzb(1,3)],elknrzb,k,m);
        A(Topology(iel,1:nknel))=A(Topology(iel,1:nknel))+h(1:nknel).';
        B(Topology(iel,1:nknel))=B(Topology(iel,1:nknel))+g(1:nknel).';
    end
    Arows(ff,:)=A;
    Brows(ff,:)=B;
    C(ff,1)=cii;
    disp(['Field point ' num2str(ff) ' of ' num2str(nfp) '   C constant: ' num2str(cii)]);
end

end