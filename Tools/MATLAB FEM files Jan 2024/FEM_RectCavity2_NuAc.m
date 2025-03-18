
% Sound field inside a rectangular cavity in 2D FEM
% Uses Matlab PDE toolbox

% The Helmholtz equation, an elliptic equation that is the time-independent
% form of the wave equation, is
%
% $-\Delta u-k^2u = 0$.
%
% Solving this equation allows us to compute the waves reflected by an
% object illuminated by incident waves that are coming from the
% left.

% Adapted from Matlab example "pdedemo2", VCH 10-2018
% See help item in Matlab: "Scattering Problem".

% 8-2019: Explicit FEM coefficient matrices

% 1-2022: The line source is implemented as a small cylinder with
% prescribed velocity. This version is therefore a radiation problem and
% there is no source term. The result seems to match better with the
% corresponding BEM calculation.


clear

% Problem Definition

% Constants
rho=1.1992;
c0=344;

% define geometry of the rectangular cavity
lx=1.2; % Heigth
ly=1; % Width

% frequency and wavenumber 
m=5;n=3;                        % Select eigenmode
fr=c0/2*sqrt((m/lx)^2+(n/ly)^2); % Use an analytical resonance
k=2*pi*fr/c0;


% line source position
p_source=[lx*0.05 ly*0.05]; % Source close to corner
% p_source=[lx/2 ly/2]; % Centered source
Rs = 0.01;  % Radius of te source as a pulsating sphere
% Normal velocity at the source. See"Fundamentals.. book by Finn J.", p.175
% and compare with the Green's function used in BEM 
V0 = besselh(1,2,k*Rs)/(4*rho*c0) ; 


% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.
c = 1; % Do not mistake for c0 (speed of sound)
a = -k^2;
f = 0; % Source term, no source (homogeneous). This is a radiation problem.

% The variables (gd,sf,ns) describe the geometry (g). They can be created, e.g.:
%
%   1-Drawing the geometry in the PDE Modeler Matlab app and exporting it to the workspace.
%   2-Creating the variables manually. The syntax is specified in the help of the "decsg" function. 
%   3-Other options: importing from CAD, etc. Check Matlab PDE toolbox help.


% Option 2: Manually define the geometry
gd =[3              1;
     4    p_source(1);
     0    p_source(2);
     lx            Rs;
     lx             0;
     0              0;
     ly             0;
     ly             0;
     0              0;
     0              0];           % Each column is a geometrical object (one rectangle here). See help of the "decsg" function.
ns = [82 82; 49 50];     % Names of the geometrical objects, characters in numeric (ASCII) form.
sf ='R1-R2';           % Represents the geometrical operation to be made.
char(ns')           % This shows the names in "ns". Must correspond to those in the operation in "sf".

g = decsg(gd,sf,ns); % Create the geometry.


% Create PDE Model
% Create a PDE Model with a single dependent variable.
numberOfPDE = 1;
model = createpde(numberOfPDE);


% Convert the geometry and append it to the pde model.
geometryFromEdges(model,g);

% Specify PDE Coefficients 
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f); % Inhomogeneous equation

% Plot the geometry and display the edge labels for use in the boundary condition definition.
figure; 
pdegplot(model,'EdgeLabels','on'); 
axis equal; grid
title 'Geometry With Edge Labels Displayed';


% Apply the boundary conditions. See help items in Matlab "Specify Boundary
% Conditions", "Solve PDEs with Nonconstant Boundary Conditions"
% Convention: exp(-jwt)
bInner = applyBoundaryCondition(model,'neumann','Edge',(1:4),'g',0,'q',0); % Hard surface
% bInner = applyBoundaryCondition(model,'neumann','Edge',(1:4),'g',0,'q',-k*1i); % This is equivalent to a (rho*c) impedance


% Source: set to prescribed velocity V0
% dp/dn=-j*w*rho*vn  (Euler)  &  dp/dn + q*p = g (PDE Matlab)  &  vn = V0 (BC) >>>>   
% g = -j*w*rho*V0 = -j*k*c*rho*V0  & q = 0
innerBCFunc = @(loc,state)-j*k*rho*c0*V0; 
bInner = applyBoundaryCondition(model,'neumann','Edge',(5:8),'g',innerBCFunc,'q',0); % Hard surface


% Create Mesh, defining mesh density
generateMesh(model,'Hmax',c0/fr/6,'Hmin',0.001); % The mesh density is determined by the last parameter. See help of "generateMesh".
figure
pdemesh(model); 
axis equal;grid
Mfem=size(model.Mesh.Nodes,2);      % Nr. of nodes
Nfem=size(model.Mesh.Elements,2);   % Nr. of elements
title(['FEM mesh: Nodes = ' num2str(Mfem) '  Elements = ' num2str(Nfem)]);
xlabel('x');ylabel('y');



% Solve using the FEM matrices:
FEM = assembleFEMatrices(model);
% p=(FEM.K+FEM.A+FEM.Q)\(FEM.F+FEM.G);  % standard form with all B.C.s and sources
% FEM.K is the stiffness matrix.
% Careful! FEM.A is the mass matrix, not M as said in 'assembleFEMatrices' help!!
% FEM.F is the source term, point source in this case
% FEM.Q contains the rho*c impedance 'Neumann' B.C. if defined, else it is zero.
% FEM.G contains the normal velocity 'Neumann' B.C. if defined, else it is zero.
% FEM.M, FEM.H, FEM.R are zero (last two correspond to Dirichlet B.C.)

% Obtain the FEM matrices for the coupling with a normal boundary velocity:
LHSfem=(FEM.K+FEM.A+FEM.Q); % Left hand side: Multiply by the unknown acoustic pressure pa
RHSfem=diag(FEM.G);         % Right hand side: Multiply by the boundary velocity (vector of ones, if already defined)
p = LHSfem\(RHSfem*ones(Mfem,1)+FEM.F);


% % Solve using PDE function (no explicit coefficient matrices):
% result = solvepde(model);
% p = result.NodalSolution;

% Plot FEM Solution
figure
pSource = i/4*besselh(0,1,k*sqrt((model.Mesh.Nodes(1,:)-p_source(1)).^2+(model.Mesh.Nodes(2,:)-p_source(2)).^2));
% pdeplot(model,'XYData',real(u),'ZData',real(p),'ColorBar','on','Mesh','off'); % Plot as a 3D surface 
pdeplot(model,'XYData',real(p),'ColorBar','on','Mesh','off');  % Plot as colormap
% pdeplot(model,'XYData',real(pSource),'ColorBar','on','Mesh','off');  % See line source only
colormap(parula)
axis equal
title(['Rectangular cavity - Frequency = ' num2str(fr) ' Hz']);



% Animate Solution to Wave Equation
% Using the solution to the Helmholtz equation, construct an animation showing
% the corresponding solution to the time-dependent wave equation.
figure
m = 30;
h = newplot; 
hf = h.Parent; 
axis tight
ax = gca;
ax.DataAspectRatio = [1 1 1];
axis off
maxu = max(abs(p));
for j = 1:m
    uu = real(exp(-j*2*pi/m*sqrt(-1))*(p));
    pdeplot(model,'XYData',uu,'ColorBar','on','Mesh','off');
    colormap(jet)
    caxis([-maxu maxu]);
    axis tight
    ax = gca;
    ax.DataAspectRatio = [1 1 1]; 
    axis off
    M(j) = getframe(hf);
end

% To play the movie 10 times, use the |movie(hf,M,10)| command.
