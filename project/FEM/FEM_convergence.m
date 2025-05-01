clear; close all; 

% Problem Definition
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c0,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
fr = 2000; 
omega = 2 * pi * fr;
k = omega/c0;
U0 = 1/(rho*c0);



c = 1;
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
%% Ground truth
mshdens = c0/fr/10;
    generateMesh(model,'Hmax',mshdens); 
   
    Mfem=size(model.Mesh.Nodes,2);      % Nr. of nodes
    Nfem=size(model.Mesh.Elements,2);   % Nr. of elements
    % VELOCITY OF PISTON
    pistBCFunc = @(loc,state)j*k*rho*c0*U0; 
    bInner = applyBoundaryCondition(model,'neumann','Edge',(EdgU0),'g',pistBCFunc,'q',0); % Hard surface
    
    % BAFFLE Surface:
    bInner = applyBoundaryCondition(model,'neumann','Edge',EdgB,'g',0,'q',0); % Hard surface
    
    
    % FREE FIELD: set to the characteristic impedance rho*c (non-reflecting).
    outerBCFunc = @(loc,state)j*k; 
    bOuter = applyBoundaryCondition(model,'neumann','Edge',EdgO,'g',0,'q',outerBCFunc); % This is equivalent to a (rho*c) impedance

 %Solve using the FEM matrices:
    % Specify PDE Coefficients 
    specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f); % Homogeneous equation

    FEM = assembleFEMatrices(model);
    % FEM.K is the stiffness matrix.
    
    LHside = (FEM.K + FEM.A+FEM.Q);
    RHside = diag(FEM.G);
    p_true = LHside\(RHside*ones(Mfem,1));
    
    Node_idx_true = model.Mesh.findNodes("region","Edge",(EdgO));
    
    x_true = model.Mesh.Nodes(1,Node_idx_true);
    y_true = model.Mesh.Nodes(2,Node_idx_true);
    [x_true,sortIdx] = sort(x_true);
    y_true = y_true(sortIdx);
    p_true = p_true(Node_idx_true(sortIdx));
    save("project/data/FEM_groundTruth.mat","fr","model","mshdens","p_true","Node_idx_true","x_true","y_true")
%% Create Mesh, defining mesh density
precis = (1:0.2:6);
p_arr = cell(1, length(precis));
N_arr = cell(1,length(precis));
X = cell(1, length(precis));
Y = cell(1, length(precis));
Mfem = zeros(1,length(precis));
for ii = 1:length(precis)
    disp("Computing:"+ii +"/"+length(precis))
    mshdens = c0/fr/precis(ii);
    generateMesh(model,'Hmax',mshdens); 
   
    Mfem(ii)=size(model.Mesh.Nodes,2);      % Nr. of nodes
    Nfem=size(model.Mesh.Elements,2);   % Nr. of elements
    
    % VELOCITY OF PISTON
    pistBCFunc = @(loc,state)j*k*rho*c0*U0; 
    bInner = applyBoundaryCondition(model,'neumann','Edge',(EdgU0),'g',pistBCFunc,'q',0); % Hard surface
    
    % BAFFLE Surface:
    bInner = applyBoundaryCondition(model,'neumann','Edge',EdgB,'g',0,'q',0); % Hard surface
    
    
    % FREE FIELD: set to the characteristic impedance rho*c (non-reflecting).
    outerBCFunc = @(loc,state)j*k; 
    bOuter = applyBoundaryCondition(model,'neumann','Edge',EdgO,'g',0,'q',outerBCFunc); % This is equivalent to a (rho*c) impedance

 %Solve using the FEM matrices:
    % Specify PDE Coefficients 
    specifyCoefficients(model,'m',0,'d',0,'c',c,'a',a,'f',f); % Homogeneous equation

    FEM = assembleFEMatrices(model);
    % FEM.K is the stiffness matrix.
    
    LHside = (FEM.K + FEM.A+FEM.Q);
    RHside = diag(FEM.G);
    p = LHside\(RHside*ones(Mfem(ii),1));
    
    Node_idx = model.Mesh.findNodes("region","Edge",(EdgO));
    x = model.Mesh.Nodes(1,Node_idx);
    y = model.Mesh.Nodes(2,Node_idx);
    [x,sortIdx] = sort(x);
    y = y(sortIdx);
    p = p(Node_idx(sortIdx));
    N_arr{ii} = Node_idx;
    X{ii} = x;
    Y{ii} = y;
    p_arr{ii} = p;
end
    save("project/data/FEM_conv_unflanged.mat","fr","model","Mfem","p_arr","Node_idx","X","Y")
%% conv error

load("FEM_groundTruth.mat")
F1 = scatteredInterpolant(x_true.',y_true.',p_true);
theta = atan2(y_,x);
theta(theta < 0) = theta(theta < 0) + 2*pi;
[theta,sortIdx] = sort(theta);

angles = [33, 55, 78, 135];
for ii= 1:length(angles)
    find()


load("FEM_groundTruthrec.mat")
F2 = scatteredInterpolant(x_true.',y_true.',p_true);
err = zeros(1, length(precis));
for ii = 1:length(err)
    true_p_interp = F(cell2mat(X(ii)),cell2mat(Y(ii)));

    err(ii) = mean(abs(true_p_interp.' - cell2mat(p_arr(ii)))./abs(true_p_interp.'));


end

figure;
loglog(Mfem,err,'-o',"LineWidth",2); 
title('Convergence of FEM (analytical is \lambda /10)' )
xlabel('Number of elements'); ylabel('Relative error')
grid

saveas(gcf, "project/figures/convergenceFEM1.svg")