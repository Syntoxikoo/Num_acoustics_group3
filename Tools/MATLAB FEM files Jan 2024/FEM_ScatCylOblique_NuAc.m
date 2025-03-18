
% Scattering by a cylinder in 2D FEM
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

clear

% Problem Definition

% Constants
rho=1.1992;
c0=344;

% frequency and wavenumber
fr=150; 
k=2*pi*fr/c0;


% define geometry of the object
Rc=1; % Radius of the cylinder
Ri=10; % Radius of the end of the domain boundary


% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.
c = 1; % Do not mistake for c0 (speed of sound)
a = -k^2;
f = 0; % Source term, no source (homogeneous). The incident waves are defined in the boundary conditions.

% The variables (gd,sf,ns) describe the geometry (g). They can be created, e.g.:
%
%   1-Drawing the geometry in the PDE Modeler Matlab app and exporting it to the workspace.
%   2-Creating the variables manually. The syntax is specified in the help
%   of the "decsg" function, and the help note "2-D Geometry Creation at Command Line"
%   3-Other options: importing from CAD, etc. Check Matlab PDE toolbox help.

% Version 1: Cylindrical outer boundary and scattering cylinder. Manually defined.
gd =[1     1;
     0     0;
     0     0;
    Ri     Rc];        % Each column is a geometrical object (circles here). See help of the "decsg" function.
ns = [69 69 ; 49 50];  % Names of the geometrical objects, characters in numeric (ASCII) form.
sf ='E1-E2';           % Represents the geometrical operation to be made.
char(ns')              % This shows the names in "ns". Must correspond to those in the operation in "sf".

% Version 2: Rectangular outer boundary and scattering cylinder. Manually defined.
% lx=Rc*30;ly=Rc*30;
% gd =[3     1;
%      4     0;
%      -lx/2 0;
%      lx/2  Rc;
%      lx/2  0;
%      -lx/2 0;
%      ly/2  0;
%      ly/2  0;
%      -ly/2 0;
%      -ly/2 0];           % Each column is a geometrical object (rectangle/circle here). See help of the "decsg" function.
% ns = [82 69; 49 49];     % Names of the geometrical objects, characters in numeric (ASCII) form.
% sf ='R1-E1';             % Represents the geometrical operation to be made.
% char(ns')                % This shows the names in "ns". Must correspond to those in the operation in "sf".


g = decsg(gd,sf,ns); % Create the geometry.


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
bOuter = applyBoundaryCondition(model,'neumann','Edge',(1:4),'g',0,'q',outerBCFunc); % This is equivalent to a (rho*c) impedance


theta=45*pi/180; % Oblique plane wave angle

% Object (cylinder) boundary: set to zero velocity (hard surface). The normal boundary velocity is set to be the opposite of the
% normal component on the surface of the prescribed incident velocity (plane wave)
% dp/dn=-j*w*rho*vn  (Euler)  &  dp/dn + q*p = g (PDE Matlab)  &  vn = -p/(rho*c)*nx = -exp(-j*k*x)/(rho*c)*nx (BC) >>>>   
% g = j*w*rho*exp(-j*k*x)/(rho*c)*nx = j*k*exp(-j*k*x)*nx  & q = 0
% innerBCFunc = @(loc,state)+j*k*exp(-1i*k*loc.x).*loc.nx; % Wave in the x-direction
innerBCFunc = @(loc,state)+j*k*exp(-1i*k*(loc.x*cos(theta)+loc.y*sin(theta))).*(loc.nx*cos(theta)+loc.ny*sin(theta)); % Version with oblique incidence
bInner = applyBoundaryCondition(model,'neumann','Edge',(5:8),'g',innerBCFunc,'q',0); % Hard surface


% Create Mesh, defining mesh density
% generateMesh(model,'Hmax',Rc/5); % The mesh density is determined by the last parameter. See help of "generateMesh".
generateMesh(model,'Hmax',c0/fr/6); % The mesh density is determined by the last parameter. See help of "generateMesh".
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
% pInc = exp(i*k*model.Mesh.Nodes(1,:)).';  % Wave in the x-direction
pInc = exp(i*k*(model.Mesh.Nodes(1,:)*cos(theta)+model.Mesh.Nodes(2,:)*sin(theta))).';   % Version with oblique incidence
subplot(1,3,1)
pdeplot(model,'XYData',real(pscat),'Mesh','off'); % Plot scattered pressure
colormap(parula);
title(['Scattered pressure around cylinder - Frequency = ' num2str(fr) ' Hz']);
subplot(1,3,2)
pdeplot(model,'XYData',real(pInc),'Mesh','off'); % Plot incident wave only
colormap(parula);
title(['Incident wave - Frequency = ' num2str(fr) ' Hz']);
subplot(1,3,3)
pdeplot(model,'XYData',real(pscat+pInc),'Mesh','off'); % Plot total pressure
colormap(parula);
title(['Total pressure around cylinder - Frequency = ' num2str(fr) ' Hz']);
% pdeplot(model,'XYData',real(pscat),'ZData',real(u),'Mesh','off'); % Plot as a 3D surface 


% Analytical solution
snodes = find(sqrt(model.Mesh.Nodes(1,:).^2+model.Mesh.Nodes(2,:).^2)<=Rc+100*eps);
ScylAxis = angle(model.Mesh.Nodes(1,snodes)+i*model.Mesh.Nodes(2,snodes))*180/pi;
[ScylAxis,iAng]=sort(ScylAxis); snodes=snodes(iAng);
% Shift analytical solution to match the incidence angle
NodesRotated=model.Mesh.Nodes(:,snodes); 
NodesRotated = [NodesRotated(1,:)*cos(-theta)-NodesRotated(2,:)*sin(-theta) ; NodesRotated(1,:)*sin(-theta)+NodesRotated(2,:)*cos(-theta)];
[ptotANA,pincANA,pscatANA]=cylscat(k,Rc,NodesRotated',150);


% Choose what to plot in the comparison
% P_FEM = pscat; P_ANA = pscatANA   % Plot scattered component only
P_FEM = pscat+pInc; P_ANA = pscatANA+pincANA;  % Plot total pressure

% plot the pressure on the surface
figure;
subplot(3,1,1)
plot(ScylAxis,real(P_FEM(snodes)),'ko',ScylAxis,real(P_ANA).','kx');grid;
title(['Scattering by a cylinder - Frequency = ' num2str(fr) ' Hz']);
xlabel('Position on the surface, deg.');ylabel('Real part of pressure (Pa)')
lax=axis; axis([-180 180 lax(3:4)]); legend('FEM','Analytical')
subplot(3,1,2)
plot(ScylAxis,abs(P_FEM(snodes)),'ko',ScylAxis,abs(P_ANA).','kx');grid;
title(['Scattering by a cylinder - Frequency = ' num2str(fr) ' Hz']);
xlabel('Position on the surface, deg.');ylabel('Pressure modulus (Pa)')
lax=axis; axis([-180 180 lax(3:4)]); legend('FEM','Analytical')
subplot(3,1,3)
plot(ScylAxis,unwrap(angle(P_FEM(snodes)))*180/pi,'ko',ScylAxis,unwrap(angle(P_ANA)).'*180/pi,'kx');grid;
title(['Scattering by a cylinder - Frequency = ' num2str(fr) ' Hz']);
xlabel('Position on the surface, deg.');ylabel('Pressure phase, deg.')
lax=axis; axis([-180 180 lax(3:4)]); legend('FEM','Analytical')




% Animate Solution to Wave Equation
% Using the solution to the Helmholtz equation, construct an animation showing
% the corresponding solution to the time-dependent wave equation.
figure
m = 10;
h = newplot; 
hf = h.Parent; 
axis tight
ax = gca;
ax.DataAspectRatio = [1 1 1];
axis off
maxu = max(abs(pscat));
for j = 1:m
    uu = real(exp(-j*2*pi/m*sqrt(-1))*(P_FEM));
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
