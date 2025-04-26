clear; close all; 

% Problem Definition


% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c0,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 

% Constants

fr = [1000 2000 4000 8000]; 
omega = 2 * pi * fr;
k = omega/c0;
% U0 = 1/(rho*c0);
U0 = 1e-3;



% * |k|, |c|, |a|, |f|: The coefficients and inhomogeneous term.

c = 1; % Do not mistake for c0 (speed of sound)
% here c corresponds to a coefficient in front of the stiffness matrix
a = -k.^2;
% here a corresponds to a coefficient in front of the mass matrix
f =0;


% The overall shape could be modified from the "FEM_model_piston.m" inside
% the geometry folder.
% save("project/FEM/geometry/flushed_piston_10meterFF.mat","gd","ns","sf"); 
load("project/FEM/geometry/flushed_pistonv2.mat")
flanged_depth = 0.1;
gd(7:8,3) =  - flanged_depth;


char(ns)


g = decsg(gd,sf,ns); % Create the geometry.


numberOfPDE = 1; 
model = createpde(numberOfPDE);


% Convert the geometry and append it to the pde model.
geometryFromEdges(model,g);


%% Plot the geometry and display the edge labels for use in the boundary condition definition.
figure(1) 
pdegplot(model,'EdgeLabels','on'); 
axis equal; grid
title 'Geometry With Edge Labels Displayed';


% Boundary condition
if flanged_depth == 0.0
    EdgU0 = 3;
    EdgB =  [1,2,4,5,6,7,8];
    EdgO = (9:12);
else
    EdgU0 = 1;
    EdgB =  (2:10);
    EdgO = (11:14);
end


%% Create Mesh, defining mesh density
p_arr = cell(1, length(fr));
X = cell(1, length(fr));
Y = cell(1, length(fr));
N_arr = cell(1,length(fr));
for ii = 1:length(fr)
    mshdens = c0/fr(ii)/6;
    generateMesh(model,'Hmax',mshdens); 
   
    Mfem=size(model.Mesh.Nodes,2);      % Nr. of nodes
    Nfem=size(model.Mesh.Elements,2);   % Nr. of elements
    if ii == 1
        figure(2)
        pdemesh(model); 
        axis equal;grid
        title(['FEM mesh: Nodes = ' num2str(Mfem) '  Elements = ' num2str(Nfem)]);
        xlabel('x');ylabel('y');
    end
    % VELOCITY OF PISTON
    pistBCFunc = @(loc,state)j*k(ii)*rho*c0*U0; 
    bInner = applyBoundaryCondition(model,'neumann','Edge',(EdgU0),'g',pistBCFunc,'q',0); % Hard surface
    
    % BAFFLE Surface:
    bInner = applyBoundaryCondition(model,'neumann','Edge',EdgB,'g',0,'q',0); % Hard surface
    
    
    % FREE FIELD: set to the characteristic impedance rho*c (non-reflecting).
    outerBCFunc = @(loc,state)j*k(ii); 
    bOuter = applyBoundaryCondition(model,'neumann','Edge',EdgO,'g',0,'q',outerBCFunc); % This is equivalent to a (rho*c) impedance

 %Solve using the FEM matrices:
    % Specify PDE Coefficients 
    specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a(ii),'f',f); % Homogeneous equation

    FEM = assembleFEMatrices(model);
    % FEM.K is the stiffness matrix.
    
    LHside = (FEM.K + FEM.A+FEM.Q);
    RHside = diag(FEM.G);
    p = LHside\(RHside*ones(Mfem,1));
    p_arr{ii} = p;
    Node_idx = model.Mesh.findNodes("region","Edge",(EdgO));
    N_arr{ii} = Node_idx;
    x = model.Mesh.Nodes(1,Node_idx);
    y = model.Mesh.Nodes(2,Node_idx);
    X{ii}= x;
    Y{ii} = y;
end
% %% Plot FEM Solution
% figure
% pdeplot(model,'XYData',20*log10(abs(p)/2e-5),'Mesh','off'); % Plot scattered pressure
% colormap(parula);
% title(['Pressure around cylinder - Frequency = ' num2str(fr) ' Hz']);

%% get directivity pattern






the_arr = cell(1,length(fr));
sor_arr = cell(1,length(fr));
figure;
for ii =1: length(fr)

    Node_idx = cell2mat(N_arr(ii));

    theta = atan2(cell2mat(Y(ii)),cell2mat(X(ii)));
    theta(theta < 0) = theta(theta < 0) + 2*pi;
    
    [theta,sortIdx] = sort(theta);
    sor_arr{ii} = sortIdx;
    the_arr{ii} = theta;
    p = cell2mat(p_arr(ii));
    spl_values = 20*log10(abs(p(Node_idx(sortIdx)))/20e-6);
    Nspl = spl_values - max(spl_values);
    polarplot(theta-pi/2, Nspl, 'LineWidth', 2, 'DisplayName', [num2str(fr(ii)) ' Hz']);

    hold on;
end

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

save("Result_FEM_flushed.mat", "the_arr","p_arr","N_arr",'fr','model', "sor_arr")