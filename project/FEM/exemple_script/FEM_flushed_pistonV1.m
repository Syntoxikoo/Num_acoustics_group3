clear; close all; 

% Problem Definition

% Constants
rho=1.21;
c0=343;
fr = 150; 
omega = 2 * pi * fr;
k = omega/c0;
U0 = 1/(rho*c0);



% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.

c = 1; % Do not mistake for c0 (speed of sound)
% here c corresponds to a coefficient in front of the stiffness matrix
a = -k^2;
% here a corresponds to a coefficient in front of the mass matrix
f =0;


% The overall shape could be modified from the "FEM_model_piston.m" inside
% the geometry folder.
% save("project/FEM/geometry/flushed_piston_10meterFF.mat","gd","ns","sf"); 
load("project/FEM/geometry/flushed_piston_10meterFF.mat")
flanged_depth = 0.0;
gd(7:8,3) =  - flanged_depth;


char(ns)


g = decsg(gd,sf,ns); % Create the geometry.


numberOfPDE = 1; 
model = createpde(numberOfPDE);


% Convert the geometry and append it to the pde model.
geometryFromEdges(model,g);

% Specify PDE Coefficients 
specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f); % Homogeneous equation

%% Plot the geometry and display the edge labels for use in the boundary condition definition.
figure(1) 
pdegplot(model,'EdgeLabels','on'); 
axis equal; grid
title 'Geometry With Edge Labels Displayed';

%% Create Mesh, defining mesh density
mshdens = c0/fr/6;
generateMesh(model,'Hmax',mshdens); 
figure(2)
pdemesh(model); 
axis equal;grid
Mfem=size(model.Mesh.Nodes,2);      % Nr. of nodes
Nfem=size(model.Mesh.Elements,2);   % Nr. of elements
title(['FEM mesh: Nodes = ' num2str(Mfem) '  Elements = ' num2str(Nfem)]);
xlabel('x');ylabel('y');
%% Boundary condition
if flanged_depth == 0.0
    EdgU0 = 3;
    EdgB =  [1,2,4,5,6,7,8];
    EdgO = (9:12);
else
    EdgU0 = 1;
    EdgB =  (2:10);
    EdgO = (10:14);
end

% VELOCITY OF PISTON
pistBCFunc = @(loc,state)j*k*rho*c0*U0; 
bInner = applyBoundaryCondition(model,'neumann','Edge',(EdgU0),'g',pistBCFunc,'q',0); % Hard surface

% BAFFLE Surface:
bInner = applyBoundaryCondition(model,'neumann','Edge',EdgB,'g',0,'q',0); % Hard surface


% FREE FIELD: set to the characteristic impedance rho*c (non-reflecting).
outerBCFunc = @(loc,state)j*k; 
bOuter = applyBoundaryCondition(model,'neumann','Edge',EdgO,'g',0,'q',outerBCFunc); % This is equivalent to a (rho*c) impedance

%% Solve using the FEM matrices:
FEM = assembleFEMatrices(model);
% FEM.K is the stiffness matrix.

LHside = (FEM.K + FEM.A+FEM.Q);
RHside = diag(FEM.G);
p = LHside\(RHside*ones(Mfem,1));


%% Plot FEM Solution
figure
pdeplot(model,'XYData',20*log10(abs(p)/2e-5),'Mesh','off'); % Plot scattered pressure
colormap(parula);
title(['Pressure around cylinder - Frequency = ' num2str(fr) ' Hz']);

%% get directivity pattern

Node_idx = model.Mesh.findNodes("region","Edge",(EdgO));


x = model.Mesh.Nodes(1,Node_idx);
y = model.Mesh.Nodes(2,Node_idx);
theta = atan2(y,x);
theta(theta < 0) = theta(theta < 0) + 2*pi;

[theta,sortIdx] = sort(theta);


figure;

spl_values = 20*log10(abs(p(Node_idx(sortIdx)))/20e-6);
Nspl = spl_values - max(spl_values);
polarplot(theta-pi/2, Nspl, 'LineWidth', 2, 'DisplayName', [num2str(fr) ' Hz']);
% 
pax = gca;
pax.ThetaZeroLocation = "top";
pax.ThetaDir = "clockwise";

grid on;
title('Directivity Pattern Comparison');
rlim([-60, 0]);
rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);
legend('Location', 'best');
hold off;

