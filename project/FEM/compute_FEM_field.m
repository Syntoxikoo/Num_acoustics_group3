clear ;close all
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c0,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 

% Constants

fr = 2000; 
omega = 2 * pi * fr;
k = omega/c0;
U0 = 1e-3 *sqrt(2);



c = 1;
a = -k.^2;
f =0;


% The overall shape could be modified from the "FEM_model_piston.m" inside
% the geometry folder.
% save("project/FEM/geometry/flushed_piston_rounded.mat","gd","ns","sf"); 
load("project/FEM/geometry/flushed_pistonv2.mat")
flanged_depth = [0.0,0.1];

p_arr = cell(1, length(flanged_depth));
X = cell(1, length(flanged_depth));
Y = cell(1, length(flanged_depth));

for ii= 1:length(flanged_depth)
    gd(7:8,3) =  - flanged_depth(ii);

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
    if flanged_depth(ii) == 0.0
        EdgU0 = 3;
        EdgB =  [1,2,4,5,6,7,8];
        EdgO = (9:12);
    else
        EdgU0 = 1;
        EdgB =  (2:10);
        EdgO = (11:14);
    end


%% Create Mesh, defining mesh density


    mshdens = c0/fr/6;
    generateMesh(model,'Hmax',mshdens); 
   
    Mfem=size(model.Mesh.Nodes,2);      % Nr. of nodes
    Nfem=size(model.Mesh.Elements,2);   % Nr. of elements

    % VELOCITY OF PISTON
    pistBCFunc = @(loc,state)1j*k*rho*c0*U0; 
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
    p = LHside\(RHside*ones(Mfem,1));
    p_arr{ii} = p;
    
    
    x = model.Mesh.Nodes(1,:);
    y = model.Mesh.Nodes(2,:);
    X{ii}= x;
    Y{ii} = y;
    if ii == 1
        
        model1 = model; 
    else
        model2 = model;
    end
end

save("project/data/Result_FEM_field.mat", "p_arr","X","Y", "model1","model2")

customColormap = [
    24/255,  23/255,  72/255;
    66/255,  65/255, 104/255;
    87/255,  85/255, 121/255;
   108/255, 107/255, 138/255;
   129/255, 129/255, 153/255;
   200/255, 200/255, 200/255;
   185/255, 125/255, 126/255;
   171/255, 102/255, 101/255;
   159/255,  75/255,  74/255;
   148/255,  49/255,  47/255;
   129/255,   1/255,   0/255
];

n_orig = size(customColormap, 1); 
n_smooth = 256; 

x_orig = linspace(0, 1, n_orig);
x_smooth = linspace(0, 1, n_smooth);

smoothColormap = [
    interp1(x_orig, customColormap(:,1), x_smooth, 'pchip');
    interp1(x_orig, customColormap(:,2), x_smooth, 'pchip');
    interp1(x_orig, customColormap(:,3), x_smooth, 'pchip')
]';
save("cmap.map","smoothColormap")

% pdeplot(model, 'XYData', real(p), 'Contour', 'on', 'ColorMap', smoothColormap)
