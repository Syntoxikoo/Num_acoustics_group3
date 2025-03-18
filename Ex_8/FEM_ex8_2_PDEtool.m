
% Modes inside a rectangular cavity in 2D FEM
% Uses Matlab PDE toolbox

% The Helmholtz equation, an elliptic equation that is the time-independent
% form of the wave equation, is
%
% $\Delta u-k^2u = 0$.
%
% Solving this equation allows us to compute the modes and egienfrequencies
% inside a cavity

% Adapted from Matlab example "pdedemo2", VCH 10-2018
% See help item in Matlab: "Scattering Problem".

% 8-2019: Exlicit FEM coefficient matrices

clear; %close all; 

% Problem Definition

% Constants
rho=1.25;
c0=343;
fr = 200; %determines the mesh density

% define geometry of the rectangular cavity
lx=10; % Heigth
ly=4;  % Width
% corresponds to the geometry in "Run_FEM_ex9_2D_Eig.m"

% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.
% Setup homogeneous eigenvalue problem:

c = 1; % Do not mistake for c0 (speed of sound)
% here c corresponds to a coefficient in front of the stiffness matrix
a = 1;
% here a corresponds to a coefficient in front of the mass matrix


% The variables (gd,sf,ns) describe the geometry (g). They can be created, e.g.:
%
%   1-Drawing the geometry in the PDE Modeler Matlab app and exporting it to the workspace.
%   2-Creating the variables manually. The syntax is specified in the help of the "decsg" function. 
%   3-Other options: importing from CAD, etc. Check Matlab PDE toolbox help.


% Option 2: Manually define the geometry
gd =[ 3 ;
      4 ;
      0 ;
     lx ;
     lx ;
      0 ;
     ly ;
     ly ;
      0 ;
      0];           % Each column is a geometrical object (one rectangle here). See help of the "decsg" function.
ns = [82 ; 49];     % Names of the geometrical objects, characters in numeric (ASCII) form.
sf ='R1';           % Represents the geometrical operation to be made.
char(ns')           % This shows the names in "ns". Must correspond to those in the operation in "sf".


g = decsg(gd,sf,ns); % Create the geometry.


% Create PDE Model
% Create a PDE Model with a single dependent variable.
numberOfPDE = 1;
model = createpde(numberOfPDE);


% Convert the geometry and append it to the pde model.
geometryFromEdges(model,g);

% Specify PDE Coefficients 
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',0); % Homogeneous equation

%% Plot the geometry and display the edge labels for use in the boundary condition definition.
figure(1) 
pdegplot(model,'EdgeLabels','on'); 
axis equal; grid
title 'Geometry With Edge Labels Displayed';

%% Create Mesh, defining mesh density
generateMesh(model,'Hmax',c0/fr/6); % The mesh density is determined by the last parameter. See help of "generateMesh".
figure(2)
pdemesh(model); 
axis equal;grid
Mfem=size(model.Mesh.Nodes,2);      % Nr. of nodes
Nfem=size(model.Mesh.Elements,2);   % Nr. of elements
title(['FEM mesh: Nodes = ' num2str(Mfem) '  Elements = ' num2str(Nfem)]);
xlabel('x');ylabel('y');

%% Solve using the FEM matrices:
FEM = assembleFEMatrices(model);
% FEM.K is the stiffness matrix.
% Careful! FEM.A is the mass matrix, not M as said in 'assembleFEMatrices' help!!

% Solve the eigenvalue problem
[X1,D]=eigs(FEM.K,FEM.A,10,'sm');

% Alternative like in "Run_FEM_ex9_2D_Eig.m"
%[X1,D]=eig(full(FEM.K),full(FEM.A));

for j=1:length(X1(1,:))
    mnorm=sqrt(X1(:,j)'*FEM.A*X1(:,j));
    X1(:,j)=X1(:,j)/mnorm;
end

[eigv,id]=sort(diag(D)); % Eigenvalue

freq = sqrt(eigv)*c0/(2*pi);

U=X1(:,id); %Eigenvector

% order mesh for patch
edof = model.Mesh.Elements';
edof = [edof(:,1),edof(:,4),edof(:,2),edof(:,5),edof(:,3),edof(:,6)];

%% Plot Modeshape
modnr=1;

figure(3);
for j=1:3 
    for i=1:3
        modnr=modnr+1;
        subplot(3,3,modnr-1)
        Ed = (U(:,modnr));
        Edabs=abs(Ed); 
        const=max(max(Edabs)); Ed=Ed/const;
        patch('faces',edof(:,1:6),...
            'vertices',[model.Mesh.Nodes',...
            zeros(length(model.Mesh.Nodes(1,:)),1)],...
            'FaceVertexCData',Ed,...
            'FaceColor','interp',...
            'LineStyle','non')
        MOD= num2str(freq(modnr));
        title(MOD)
        axis off
        axis equal
    end 
end
 