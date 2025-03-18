
% Sound field inside a cylindrical cavity in axisimetrical FEM, based on 2D
% Uses Matlab PDE toolbox
% Compare with "TestRadCylinder.m" in axisymmetrical BEM 

% The Helmholtz equation, an elliptic equation that is the time-independent
% form of the wave equation, is
%
% $-\Delta u-k^2u = 0$.
%
% Solving this equation allows us to compute the waves reflected by an
% object illuminated by incident waves that are coming from the
% left.

% Adapted from Matlab example "radioactiveRod.m", VCH 10-2018
% See help item in Matlab: "Heat Distribution in Circular Cylindrical Rod": https://se.mathworks.com/help/pde/examples/heat-distribution-in-a-circular-cylindrical-rod.html

% 8-2019: Explicit FEM coefficient matrices


clear

% Problem Definition

% Constants
rho=1.1992;
c0=344;

% define geometry of the rectangular cavity
Hc=1; % Heigth
Rc=0.4;   % Radius

% frequency and wavenumber 
% Reference: "Vibration and Sound", Morse, p.398
Alpha0n = [0 1.2197 2.2331 3.2383 4.2411]; % Only radial and axial waves, all axisymmetrical
nz=2; Alpha=Alpha0n(4);            % Select eigenmode
fr=c0/2*sqrt((nz/Hc)^2+(Alpha/Rc)^2); % Use an analytical resonance
k=2*pi*fr/c0;


% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.
% c = 1; % Do not mistake for c0 (speed of sound)
% a = -k^2;
c = @(region,state)+region.x; % Do not mistake for c0 (speed of sound)
a = @(region,state)-k^2*region.x;
f = 0; % Source term, no source (homogeneous). The incident waves are defined in the boundary conditions.


% The variables (gd,sf,ns) describe the geometry (g). They can be created, e.g.:
%
%   1-Drawing the geometry in the PDE Modeler Matlab app and exporting it to the workspace.
%   2-Creating the variables manually. The syntax is specified in the help of the "decsg" function. 
%   3-Other options: importing from CAD, etc. Check Matlab PDE toolbox help.


% Option 2: Manually define the geometry
gd =[3      3      3;
     4      4      4;
     0      0      0;
     Rc     Rc     0.8*Rc;
     Rc     Rc     0.8*Rc;
     0      0      0;
     
     Hc    0.8*Hc  Hc;  
     Hc    0.8*Hc  Hc;
     0     0.2*Hc  0;
     0     0.2*Hc  0];           % Each column is a geometrical object (one rectangle here). See help of the "decsg" function.
ns = [82 82 82; 49 50 51];     % Names of the geometrical objects, characters in numeric (ASCII) form.
sf ='R1+R2+R3';           % Represents the geometrical operation to be made.
char(ns')           % This shows the names in "ns". Must correspond to those in the operation in "sf".


[g0,bt] = decsg(gd,sf,ns); % Create the geometry. 
[g,bt2] = csgdel(g0,bt);   % Remove internal boundaries, if any


% Create PDE Model
% Create a PDE Model with a single dependent variable.
numberOfPDE = 1;
model = createpde(numberOfPDE);


% Convert the geometry and append it to the pde model.
geometryFromEdges(model,g);


% Specify PDE Coefficients 
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f);


% Plot the geometry and display the edge labels for use in the boundary condition definition.
figure; 
pdegplot(model,'EdgeLabels','on','SubdomainLabels','on')
axis equal; grid
title 'Geometry With Edge Labels Displayed';


% Apply the boundary conditions. See help items in Matlab "Specify Boundary
% Conditions", "Solve PDEs with Nonconstant Boundary Conditions"
% Convention: exp(-jwt)
% Define a hard surface on part of the boundary:
bWall = applyBoundaryCondition(model,'neumann','Edge',([1 2 3 9 10]),'g',0,'q',0);

% Define an impedance on part of the boundary (rho*c is non-reflecting):
% dp/dn=-j*w*rho*vn  (Euler)  &  c*dp/dn + q*p = g (PDE Matlab)  &  p/vn = -(rho*c) (BC) (normal opposite to propagation direction) >>>>   
% g = 0  & q = -dp/dn*(1/p) = j*w*rho*vn*(1/p) = -j*w*rho/(rho*c) = -j*k*c  (c=x)
ZWallBCFunc = @(loc,state)-1i*k.*loc.x; 
ZWall = applyBoundaryCondition(model,'neumann','Edge',([4 5]),'g',0,'q',ZWallBCFunc);

% Define a normal velocity on part of the boundary:
% dp/dn=-j*w*rho*vn  (Euler)  &  c*dp/dn + q*p = g (PDE Matlab)  &  vn = 1/(rho*c0) (BC) >>>>   
% g = -j*w*rho/(rho*c0)*c (c=x) >>  g = -j*k*c (c=x)  &   q = 0
RadBCFunc = @(loc,state)-1i*k*loc.x.*loc.ny; % vn=1/(rho*c)  - The "c" constant from the Diff. Eq. is also present in the Neumann BC. See "Specify Boundary Conditions"
% RadBCFunc = @(loc,state)-1i*k*c0*rho*loc.x; % vn=1   - The "c" constant from the Diff. Eq. is also present in the Neumann BC. See "Specify Boundary Conditions"
bRad = applyBoundaryCondition(model,'neumann','Edge',([6 7]),'g',RadBCFunc,'q',0);



% Create Mesh, defining mesh density
generateMesh(model,'Hmax',c0/fr/6); % The mesh density is determined by the last parameter. See help of "generateMesh".
figure
pdemesh(model); 
axis equal;grid
Mfem=size(model.Mesh.Nodes,2);      % Nr. of nodes
Nfem=size(model.Mesh.Elements,2);   % Nr. of elements
title(['FEM mesh: Nodes = ' num2str(Mfem) '  Elements = ' num2str(Nfem)]);
xlabel('rho');ylabel('z');



% Solve using the FEM matrices:
FEM = assembleFEMatrices(model);
% p=(FEM.K+FEM.A+FEM.Q)\(FEM.F+FEM.G);  % standard form with all B.C.s and sources
% FEM.K is the stiffness matrix.
% Careful! FEM.A is the mass matrix, not M as said in 'assembleFEMatrices' help!!
% FEM.F is the source term, zero in this case
% FEM.Q contains the rho*c impedance 'Neumann' B.C. if defined, else it is zero.
% FEM.G contains the normal velocity 'Neumann' B.C. if defined, else it is zero.
% FEM.M, FEM.H, FEM.R are zero (last two correspond to Dirichlet B.C.)

% Obtain the FEM matrices for the coupling with a normal boundary velocity:
LHSfem=(FEM.K+FEM.A+FEM.Q); % Left hand side: Multiply by the unknown acoustic pressure pa
RHSfem=diag(FEM.G);         % Right hand side: Multiply by the boundary velocity (vector of ones, if already defined)
p = LHSfem\(RHSfem*ones(Mfem,1));


% Solve using PDE function (no explicit coefficient matrices):
% result = solvepde(model);
% p = result.NodalSolution;

% Plot FEM Solution
figure
pdeplot(model,'XYData',abs(p),'ZData',abs(p),'ColorBar','on','Mesh','on'); % Plot as a 3D surface 
% pdeplot(model,'XYData',real(p),'ColorBar','on','Mesh','off');  % Plot as colormap
colormap(parula);grid
title(['Cylinder (FEM), ka = ' num2str(k*Rc) ' , f = ' num2str(fr) ' Hz'])
xlabel('r (m)');ylabel('z (m)')


% Animate Solution to Wave Equation
% Using the solution to the Helmholtz equation, construct an animation showing
% the corresponding solution to the time-dependent wave equation.
figure
m = 60;
h = newplot; 
hf = h.Parent; 
axis tight
ax = gca;
ax.DataAspectRatio = [1 1 1];
axis off
maxu = max(abs(p));
for j = 1:m
    uu = real(exp(-j*2*pi/m*sqrt(-1))*(p));
%     pdeplot(model,'XYData',uu,'ColorBar','on','Mesh','off');
    pdeplot(model,'XYData',uu,'ZData',uu,'ColorBar','on','Mesh','on'); % Plot as a 3D surface 
    colormap(parula);grid
    lax=axis; axis([lax(1:4) -maxu maxu]);
%     caxis([-maxu maxu]);
%     axis tight
    ax = gca;
%     ax.DataAspectRatio = [1 1 1]; 
    title(['Rectangular cavity - k*a = ' num2str(k*Rc) ' Hz']);
%     axis off
    M(j) = getframe(hf);
end

% To play the movie 10 times, use the |movie(hf,M,10)| command.
