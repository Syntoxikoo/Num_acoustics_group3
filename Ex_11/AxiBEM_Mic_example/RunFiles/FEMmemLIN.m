function [Am,Bm,rhs]=FEMmemLIN(rmem)

% [Am,Bm,rhs]=FEMmemLIN(rmem);
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
%     - rmem : vector, node coordinates (radii), from the center to the
%              rim. Linear one-dimensional elements are used.
%
% Outputs:
%     - Am,Bm: coefficient matrices.
%     - rhs  : right hand side matrix, excitation.
%
% Reference: "The Finite Element Method using MATLAB", Y.W.Kwon, H.Bang

% The weak Galerkin formulation of the FEM is used.
% The membrane equation is as in Morse (Vibration and Sound, Ch:5, Sec 17), Robey1954, Zuckerwar 1978
% If the laplacian is expressed in polar coordinates. The 1/r*d(disp)/dr
% term vanishes when integrating by parts. The surface differential is
% 2*pi*r*dr.
% Vicente Cutanda 07-2003

ng=50; % number of Gauss points
[bp,wf]=gaussrule(ng);

nmem=length(rmem);
Am=zeros(nmem);
Bm=zeros(nmem);
for ii=1:nmem-1 % loop over every element
   % coordinates of the nodes
   ri1=rmem(ii);ri2=rmem(ii+1); 
   rr=(ri1+ri2)/2+(ri2-ri1)/2*bp;
   weight=(ri2-ri1)/2*wf;
   
   % Shape functions and their derivatives, no local coordinates used
   H=[(ri2-rr)./(ri2-ri1) (rr-ri1)./(ri2-ri1)];
   dH=[-1./(ri2-ri1) 1./(ri2-ri1)];
   
   % building of the element FEM matrices (2x2) - weak formulation
   for iif=1:2
      for iic=1:2
         AMtmp=-(dH(:,iif).*dH(:,iic))*2*pi.*rr;
         BMtmp=(H(:,iif).*H(:,iic))*2*pi.*rr;
         
         AM(iif,iic)=weight'*AMtmp;
         BM(iif,iic)=weight'*BMtmp;
       end
   end
   
   % add element FEM matrices to the global FEM matrices
   Am(ii:ii+1,ii:ii+1)=Am(ii:ii+1,ii:ii+1)+AM;
   Bm(ii:ii+1,ii:ii+1)=Bm(ii:ii+1,ii:ii+1)+BM;
end

rhs=Bm; % the right hand side matrix is equal to Bm
