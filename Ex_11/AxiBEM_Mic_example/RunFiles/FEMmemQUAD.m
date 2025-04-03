function [Am,Bm,rhs,rnode,rtopo]=FEMmemQUAD(rmem)

% [Am,Bm,rhs,rnode,rtopo]=FEMmemQUAD(rmem);
% 
% 
% Calculates the FEM coefficient matrices for a circular membrane
% using Galerkin formulation. No boudary conditions are applied.
%
% The resulting equation is: Tm*(Am+K2*Bm)*[disp]=rhs*[pressure]
% where [disp] and [pressure] are the displacement and pressure 
% column vectors. K2 is the square of the wavenumber on the membrane (1/m)
% and Tm is the membrane tension (N/m). Boundary conditions must be set.
% 
% The formulation is fully axisymmetrical.
% 
% Inputs:
%     - rmem : vector, coordinates (radii) of the elements' connecting nodes,
%              from the center to the rim. Quadratic one-dimensional elements
%              are used. The mid node of each element is added in the output.
%
% Outputs:
%     - Am,Bm: coefficient matrices.
%     - rhs  : right hand side matrix, excitation.
%     - rnode: node radii, with the added midpoints.
%     - rtopo: topology matrix, relating nodes and elements.
%
% Reference: "The Finite Element Method using MATLAB", Y.W.Kwon, H.Bang

% Vicente Cutanda 07-2003

% The weak Galerkin formulation of the FEM is used.
% The membrane equation is as in Morse (Vibration and Sound, Ch:5, Sec 17), Robey1954, Zuckerwar 1978
% If the laplacian is expressed in polar coordinates. The 1/r*d(disp)/dr
% term vanishes when integrating by parts. The surface differential is
% 2*pi*r*dr.

ng=50; % number of Gauss points
[bp,wf]=gaussrule(ng);

nmem=length(rmem)-1; % number of elements
Am=zeros(2*nmem+1);
Bm=zeros(2*nmem+1);
rnode=zeros(2*nmem+1,1);
rtopo=zeros(nmem,3);
for ii=1:nmem % loop over every element
   % coordinates of the nodes
   ri1=rmem(ii);ri2=(rmem(ii)+rmem(ii+1))/2;ri3=rmem(ii+1);
   rnode(ii*2-1:ii*2+1,1)=[ri1;ri2;ri3];
   rtopo(ii,:)=ii*2-1:ii*2+1;
   
   rr=(ri1+ri3)/2+(ri3-ri1)/2*bp;
   weight=(ri3-ri1)/2*wf;
   
   % Shape functions and their derivatives, no local coordinates used
   H=2/(ri3-ri1)^2*[(rr-ri2).*(rr-ri3) -2*(rr-ri1).*(rr-ri3) (rr-ri1).*(rr-ri2)];
   dH=2/(ri3-ri1)^2*[2*rr-ri2-ri3 -2*(2*rr-ri1-ri3) 2*rr-ri1-ri2];
       
   % building of the element FEM matrices (2x2) - weak formulation
   for iif=1:3
      for iic=1:3
         AMtmp=-(dH(:,iif).*dH(:,iic))*2*pi.*rr;
         BMtmp=(H(:,iif).*H(:,iic))*2*pi.*rr;
         
         AM(iif,iic)=weight'*AMtmp;
         BM(iif,iic)=weight'*BMtmp;
       end
   end
   
   % add element FEM matrices to the global FEM matrices
   Am(ii*2-1:ii*2+1,ii*2-1:ii*2+1)=Am(ii*2-1:ii*2+1,ii*2-1:ii*2+1)+AM;
   Bm(ii*2-1:ii*2+1,ii*2-1:ii*2+1)=Bm(ii*2-1:ii*2+1,ii*2-1:ii*2+1)+BM;
end

rhs=Bm; % the right hand side matrix is equal to Bm
