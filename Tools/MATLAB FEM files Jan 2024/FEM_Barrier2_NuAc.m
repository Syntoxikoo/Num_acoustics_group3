
% Scattering by a barrier on a plane in 2D FEM
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
% See help item in Matlab: "Scattering Problem": https://www.mathworks.com/help/pde/ug/scattering-problem.html

% 8-2019: Explicit FEM coefficient matrices

% 1-2022: The line source is implemented as a small cylinder with
% prescribed velocity. This version is therefore a radiation problem and
% there is no source term. The result seems to match better with the
% corresponding BEM calculation. The source definition has been validated
% against an analytical expression in another calculation.

clear

% Problem Definition

% Constants
rho=1.1992;
c0=344;

% frequency and wavenumber
fr=150; 
k=2*pi*fr/c0;


% define geometry of the object
Bhgt=3;     % barrier height
Bwdt=0.2;   % barrier width
Ri=15;      % Radius of the end of the domain boundary

% line source position
% p_source=[-Bwdt*10 Bhgt/5]; 
p_source=[-10 0.5];
Rs = 0.1;  % Radius of te source as a pulsating cylinder
% Normal velocity at the source. See "Fundamentals..." book by Finn J.", p.175
% and compare with the Green's function used in BEM 
V0 = besselh(1,2,k*Rs)/(4*rho*c0); 


% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.
c = 1; % Do not mistake for c0 (speed of sound)
a = -k^2;
f = 0; % Source term, no source (homogeneous). This is a radiation problem.

% The variables (gd,sf,ns) describe the geometry (g). They can be created, e.g.:
%
%   1-Drawing the geometry in the PDE Modeler Matlab app and exporting it to the workspace.
%   2-Creating the variables manually. The syntax is specified in the help of the "decsg" function. 
%   3-Other options: importing from CAD, etc. Check Matlab PDE toolbox help.

% Rectangular outer boundary and scattering cylinder. Manually defined.
gd =[3     1      3        1;
     4     0      4        p_source(1);
     -Ri*2 0    -Bwdt/2    p_source(2);
     Ri*2  Ri    Bwdt/2    Rs;
     Ri*2  0     Bwdt/2    0;
     -Ri*2 0    -Bwdt/2    0;
     Ri*2  0     Bhgt      0;
     Ri*2  0     Bhgt      0;
     0     0      0        0;
     0     0      0        0];    % Each column is a geometrical object (rectangles/circles here). See help of the "decsg" function.
ns = [82 69 69 82; 49 49 50 50];%[82 69 69; 49 49 50];     % Names of the geometrical objects, characters in numeric (ASCII) form.
sf ='E1*R1-E2-R2';             % Represents the geometrical operation to be made.
char(ns')                % This shows the names in "ns". Must correspond to those in the operation in "sf".


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

% Outer boundary: set to the characteristic impedance rho*c (non-reflecting).
% dp/dn=-j*w*rho*vn  (Euler)  &  dp/dn + q*p = g (PDE Matlab)  &  p/vn = -(rho*c) (BC) (normal opposite to propagation direction) >>>>   
% g = 0  & q = -dp/dn*(1/p) = j*w*rho*vn*(1/p) = -j*w*rho/(rho*c) = -j*k
outerBCFunc = @(loc,state)-j*k; 
bOuter = applyBoundaryCondition(model,'neumann','Edge',(6:8),'g',0,'q',outerBCFunc); % This is equivalent to a (rho*c) impedance

% Object (barrier+ground) boundary: set to zero velocity (hard surface). 
bInner = applyBoundaryCondition(model,'neumann','Edge',(1:5),'g',0,'q',0); % Hard surface


% Source: set to prescribed velocity V0
% dp/dn=-j*w*rho*vn  (Euler)  &  dp/dn + q*p = g (PDE Matlab)  &  vn = V0 (BC) >>>>   
% g = -j*w*rho*V0 = -j*k*c*rho*V0  & q = 0
innerBCFunc = @(loc,state)-j*k*rho*c0*V0; 
bInner = applyBoundaryCondition(model,'neumann','Edge',(9:12),'g',innerBCFunc,'q',0); % Hard surface



% Create Mesh, defining mesh density
generateMesh(model,'Hmax',c0/fr/6,'Hmin',0.01); % The mesh density is determined by the last parameter. See help of "generateMesh".
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


% Solve using PDE function (no explicit coefficient matrices):
% result = solvepde(model);
% p = result.NodalSolution; % scattered pressure

% Plot FEM Solution
figure
pSource = i/4*besselh(0,1,k*sqrt((model.Mesh.Nodes(1,:)-p_source(1)).^2+(model.Mesh.Nodes(2,:)-p_source(2)).^2)) + ...
          i/4*besselh(0,1,k*sqrt((model.Mesh.Nodes(1,:)-p_source(1)).^2+(model.Mesh.Nodes(2,:)+p_source(2)).^2));  % Include mirror source
p_dB=20*log10(abs(p/20e-6/sqrt(2)));
pRel_dB=20*log10(abs(p./pSource));
subplot(2,1,1)
pdeplot(model,'XYData',p_dB,'ColorBar','on','Mesh','off'); 
xlabel('Distance (m)'); ylabel('Height (m)');
axis equal; colormap(parula)
title(['SPL (dB) - Frequency = ' num2str(fr) ' Hz']);
subplot(2,1,2)
pdeplot(model,'XYData',pRel_dB,'ColorBar','on','Mesh','off');
xlabel('Distance (m)'); ylabel('Height (m)');
axis equal; colormap(parula)
title(['SPL rel. free field (dB) - Frequency = ' num2str(fr) ' Hz']);



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


