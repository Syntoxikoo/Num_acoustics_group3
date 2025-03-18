function [A,B,C]=BEMEquat0(rzb,Topology,k,m,varargin)

%  Calculate C constants and A,B matrices for Helmholtz equation:
%
%  [A,B,C]=BEMEquat0(rzb,Topology,k,m,{chiefpoints,planeON})
%
%  Axisymetric BEM calculation.
%  
%  Input parameters (in SI units):
%  
%    -rzb     : Geometry matrix
%               one row for each node
%               rho-coordinate in column 1
%               z-coordinate in column 2
%               body number in column 3
%    -Topology: Topology matrix
%               One row for each element
%               Each row contains the global node number
%               in the column no. corresponding to the local
%               node number. Last column contains the body number.
%               (If left empty it is generated automatically)
%    -k       : Wave number (k=2*pi*f/c) (real scalar)
%    -m       : Circumferential mode numbers (vector of integers>=0)
%    -chiefpoints: Like 'rzb', but contains CHIEF points instead
%               One row for each chief point. Include a body number 1.
%    -planeON : if 1, it assumes a reflecting plane at z=0, using Green's function
%
%  Output parameters:
%
%    -A :       A coefficient matrix, the indexes are:
%               A(M nodes,M nodes,m terms)
%    -B :       B coefficient matrix, the indexes are:
%               B(M nodes,M nodes,m terms)
%    -C :       real column vector containing the C-constants
%
% When non-continuous boundary conditions exist, the columns in
% matrix B can be expanded. An imaginary part must be added to
% the body numbers in 'topology' to indicate what boundary area
% the element belongs to (1,2,3...). The function 'bound2D' can be
% used to set up the equations and define the boundary areas.(VC)

% VCH 5-2016: A reflecting plane at z=0 can be defined.
% Could be improved by implementing the new Green's function into the
% integral functions, rather than calling them twice. An impedance term
% could also be included, like in 2DBEM.

% VCH 11-2018: B-matrix expansion, BC splitting. 

M=size(rzb,1);
[nel, ncols]=size(Topology);
nknel=ncols-1;
if nargin<5
   chiefpoints=[];
else
   chiefpoints=varargin{1};
end
nchiefp=size(chiefpoints,1);

if nargin>5 && varargin{2}==1
   operate='y';
   if any(rzb(:,2)<0)
       error('Objects cannot reach under z=0 when an impedance plane is defined');
   end
else
   operate='n';
end

[rzbnd,Topologynd,rzdum,topodum,rznodum]=nodummy(rzb,Topology,operate);
lndu=size(rzbnd,1);

C=zeros(lndu+nchiefp,1);
A=zeros(lndu+nchiefp,lndu,length(m));
B=zeros(lndu+nchiefp,lndu,length(m));



% Prepare expansion of columns in matrix B, if the boundary conditions are not continuous.
% CConst=2*pi*(1+sign(xybnd(:,3)))/2; % 2*pi or 0 for exterior/interior domain.
% A=zeros(lndu+nchiefp,lndu);
if isreal(Topology) | operate=='y' % calculations with impedance plane are not combined with B matrix extension
    BCtopo=Topologynd;
    B=A;
else 
    BCelem=imag(Topology(:,end));
    Topologynd=real(Topology);
    Topology=real(Topology);
    % see help in function 'bound' (VC)
    [BCtopo,BCnodeA,BCnodeB]=bound2D(rzb,Topology,BCelem,'n');
    B=zeros(lndu+nchiefp,size(BCnodeB,1),length(m));
end



% Find 'A' and 'B' matrices
%for ii=1:M+nchiefp
for ii=1:lndu+nchiefp
  if ii <= lndu
     przb=rzbnd(ii,:);
  else
     przb=chiefpoints(ii-lndu,:);
  end
  cii=4*pi*(1+sign(przb(3)))/2; % 4*pi or 0 for exteror/interior domain.(VC)
  Am=zeros(lndu,length(m));
  if exist('BCnodeB')
      Bm=zeros(size(BCnodeB,1),length(m));
  else
      Bm=zeros(lndu,length(m));
  end
  for iel=1:nel
      elknrzb(1:nknel,:)=rzb(Topology(iel,1:nknel),:);
      % Singular part
      [g,h,cjj] = intF2(przb,elknrzb,m);
      if ii <= lndu
          cii=cii+cjj*(1+rzdum(ii));
      end
      if topodum(iel)~=0
          Am(Topologynd(topodum(iel),1:nknel),:)=Am(Topologynd(topodum(iel),1:nknel),:)+h(1:nknel,:);
          Bm(BCtopo(topodum(iel),1:nknel),:)=Bm(BCtopo(topodum(iel),1:nknel),:)+g(1:nknel,:);
      end
      if operate=='y'
          [g,h] = intF2([przb(1) -przb(2) 1],elknrzb,m);
          if topodum(iel)~=0
              Am(Topologynd(topodum(iel),1:nknel),:)=Am(Topologynd(topodum(iel),1:nknel),:)+h(1:nknel,:);
              Bm(BCtopo(topodum(iel),1:nknel),:)=Bm(BCtopo(topodum(iel),1:nknel),:)+g(1:nknel,:);
          end
      end
      
      % Oscillating part
      [g,h] = intF1(przb,elknrzb,k,m);
      if topodum(iel)~=0
          Am(Topologynd(topodum(iel),1:nknel),:)=Am(Topologynd(topodum(iel),1:nknel),:)+h(1:nknel,:);
          Bm(BCtopo(topodum(iel),1:nknel),:)=Bm(BCtopo(topodum(iel),1:nknel),:)+g(1:nknel,:);
      end
      if operate=='y'
          [g,h] = intF1([przb(1) -przb(2) 1],elknrzb,k,m);
          if topodum(iel)~=0
              Am(Topologynd(topodum(iel),1:nknel),:)=Am(Topologynd(topodum(iel),1:nknel),:)+h(1:nknel,:);
              Bm(BCtopo(topodum(iel),1:nknel),:)=Bm(BCtopo(topodum(iel),1:nknel),:)+g(1:nknel,:);
          end
      end
  end
  if ii/100==round(ii/10)
%       disp(sprintf('Row %g/%g C constant = %1.15e',ii,M+nchiefp,cii/4/pi));
      disp([' AxiBEM calculation, row ' num2str(ii) ' of ' num2str(lndu+nchiefp) ', C constant = ' num2str(cii/pi) ' * pi'])
  end
  
  if ii <= lndu
      Am(ii,:)=Am(ii,:)-cii;
      C(ii)=cii;
  end
  A(ii,:,:)=Am(:,:);
  B(ii,:,:)=Bm(:,:);
end

end
