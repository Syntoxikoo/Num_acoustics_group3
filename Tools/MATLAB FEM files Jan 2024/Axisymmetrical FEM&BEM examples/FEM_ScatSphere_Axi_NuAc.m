
% Scattering by a sphere in axisymmetrical FEM
% Uses Matlab PDE toolbox with a 2D definition

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

% Wavenumber and frequency
k=5;
fr=k*c0/(2*pi); 


% define geometry of the object
Rs=1;  % Radius of the sphere
Ri=10; % Radius of the end of the domain boundary


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

% Spherical outer boundary and scattering sphere. Manually defined.
gd =[3     1    1;
     4     0    0;
     -Ri*2 0    0;
     0     Rs   Ri;
     0     0    0;
     -Ri*2 0    0;
     Ri*2  0    0;
     Ri*2  0    0;
     -Ri*2 0    0;
     -Ri*2 0    0];           % Each column is a geometrical object (rectangle/circle here). See help of the "decsg" function.
ns = [82 69 69; 49 49 50];     % Names of the geometrical objects, characters in numeric (ASCII) form.
sf ='E2-E1-R1';             % Represents the geometrical operation to be made.
char(ns')                % This shows the names in "ns". Must correspond to those in the operation in "sf".


% g = decsg(gd,sf,ns); % Create the geometry.
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
bWall = applyBoundaryCondition(model,'neumann','Edge',([1 2]),'g',0,'q',0);

% Outer boundary: set to the characteristic impedance rho*c (non-reflecting).
% dp/dn=-j*w*rho*vn  (Euler)  &  c*dp/dn + q*p = g (PDE Matlab)  &  p/vn = -(rho*c) (BC) (normal opposite to propagation direction) >>>>   
% g = 0  & q = -dp/dn*(1/p) = j*w*rho*vn*(1/p) = -j*w*rho/(rho*c) = -j*k*c  (c=x)
outerBCFunc = @(loc,state)-1i*k.*loc.x; 
bOuter = applyBoundaryCondition(model,'neumann','Edge',([5 6]),'g',0,'q',outerBCFunc); % This is equivalent to a (rho*c) impedance


% Object (sphere) boundary: set to zero velocity (hard surface). The normal boundary velocity is set to be the opposite of the
% normal component on the surface of the prescribed incident velocity (plane wave)
% dp/dn=-j*w*rho*vn  (Euler)  &  c*dp/dn + q*p = g (PDE Matlab)  &  vn = -p/(rho*c)*nx = -exp(-j*k*x)/(rho*c)*nx (BC) >>>>   
% g = j*w*rho*exp(-j*k*x)/(rho*c)*nx*c = j*k*exp(-j*k*x)*nx*c   (c=x)   &     q = 0
innerBCFunc = @(loc,state)+1i*k*exp(-1i*k*loc.y).*loc.ny.*loc.x;   % The "c" constant from the Diff. Eq. is also present in the Neumann BC. See "Specify Boundary Conditions"
bInner = applyBoundaryCondition(model,'neumann','Edge',([3 4]),'g',innerBCFunc,'q',0); % Hard surface


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
pscat = LHSfem\(RHSfem*ones(Mfem,1));


% Solve using PDE function (no explicit coefficient matrices):
% result = solvepde(model);
% pscat = result.NodalSolution; % scattered pressure

% Plot FEM Solution
figure
pInc = exp(i*k*model.Mesh.Nodes(2,:)).';
subplot(1,3,1)
pdeplot(model,'XYData',real(pscat),'Mesh','off'); % Plot scattered pressure
colormap(parula);
title(['Scattered pressure around sphere - Frequency = ' num2str(fr) ' Hz']);
subplot(1,3,2)
pdeplot(model,'XYData',real(pInc),'Mesh','off'); % Plot incident wave only
colormap(parula);
title(['Incident wave - Frequency = ' num2str(fr) ' Hz']);
subplot(1,3,3)
pdeplot(model,'XYData',real(pscat+pInc),'Mesh','off'); % Plot total pressure
colormap(parula);
title(['Total pressure around sphere - Frequency = ' num2str(fr) ' Hz']);
% pdeplot(model,'XYData',real(pscat),'ZData',real(u),'Mesh','off'); % Plot as a 3D surface 


% Analytical solution on the boundary
ScylAxis_AN = linspace(0,180,150);
pAN=zeros(length(ScylAxis_AN),1);
for ll=1:length(ScylAxis_AN)
    pAN(ll) = PlaneWaveScatSphere(k,Rs,Rs,ScylAxis_AN(ll)*pi/180,1e-6); 
end
% pAN=conj(pAN);
pIncAN = conj(exp(i*k*Rs*cos(ScylAxis_AN*pi/180))).';


% plot the pressure on the surface
% Choose what to plot in the comparison
P_FEM = pscat + pInc; P_ANA = pAN;  textA= '(scattered+incident)';   % Plot total pressure
% P_FEM = pscat; P_ANA = pAN - pIncAN;  textA= '(scattered)';   % Plot scattered component only
snodes = find(abs(model.Mesh.Nodes(1,:)+1i*model.Mesh.Nodes(2,:))<=Rs+100*eps);
ScylAxis = angle(model.Mesh.Nodes(1,snodes)+1i*model.Mesh.Nodes(2,snodes))*180/pi+90;
figure;
subplot(3,1,1)
plot(ScylAxis,real(P_FEM(snodes)),'ko',ScylAxis_AN,real(P_ANA).','k.');grid;
title(['Pressure ' textA ' around a sphere - f = ' num2str(fr) ' Hz, ka = ' num2str(k*Rs)]);
xlabel('Position on the surface, deg.');ylabel('Real part of pressure (Pa)')
lax=axis; axis([0 180 lax(3:4)]); legend('FEM','Analytical')
subplot(3,1,2)
plot(ScylAxis,abs(P_FEM(snodes)),'ko',ScylAxis_AN,abs(P_ANA).','k.');grid;
xlabel('Position on the surface, deg.');ylabel('Pressure modulus (Pa)')
lax=axis; axis([0 180 lax(3:4)]); legend('FEM','Analytical')
subplot(3,1,3)
% plot(ScylAxis,unwrap(angle(P_FEM(snodes)))*180/pi,'ko',ScylAxis_AN,unwrap(angle(P_ANA)).'*180/pi,'kx');grid;
plot(ScylAxis,(angle(P_FEM(snodes)))*180/pi,'ko',ScylAxis_AN,(angle(P_ANA)).'*180/pi,'k.');grid;
xlabel('Position on the surface, deg.');ylabel('Pressure phase, deg.')
lax=axis; axis([0 180 lax(3:4)]); legend('FEM','Analytical')




% Animate Solution to Wave Equation
% Using the solution to the Helmholtz equation, construct an animation showing
% the corresponding solution to the time-dependent wave equation.
figure
m = 50;
h = newplot; 
hf = h.Parent; 
axis tight
ax = gca;
ax.DataAspectRatio = [1 1 1];
axis off
maxu = max(abs(pscat));
for j = 1:m
    uu = real(exp(j*2*pi/m*sqrt(-1))*(P_FEM));
    pdeplot(model,'XYData',uu,'ColorBar','off','Mesh','off');
    colormap(jet)
    caxis([-maxu maxu]);
    axis tight
    ax = gca;
    ax.DataAspectRatio = [1 1 1]; 
    axis off
    M(j) = getframe(hf);
end

% To play the movie 10 times, use the |movie(hf,M,10)| command.
