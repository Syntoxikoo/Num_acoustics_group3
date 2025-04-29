% 2D BEM influence of mesh and field point density
% -------------------- INPUT DATA ---------------------
 
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 
% General parameters:

fr =  2000;
vampl=1;             % Amplitude of the diaphragm movement (m/s) 
betaP = NaN;           % normalised admittance of the plane, at k.
% Field points parameters
thetamax=2*pi;        % Maximum angle for representation
Rmed=1;               % Radius of the arc of field points
beta=0; 

% DEFINITION OF DOMAIN GEOMETRY AND EXCITATION
piston_rad = 0.1;
pist_depth = 0.1;
piston_ctr = 0;
baffle_rad = 0.3;
baffle_z = 0;
thickness = 0.3; % changing thickness seems to have no effect since interior pb is omitted





% -------------- Low mesh density -------------

el_wl=1 * fr/c;
espac=1/el_wl;       
fp_spc_deg = asin(espac/Rmed)*180/pi;      
kp=2*pi*fr/c;        


segments=[-baffle_rad baffle_z -piston_rad baffle_z 5 0 el_wl;
-piston_rad baffle_z -piston_rad -pist_depth 5 0 el_wl;
          -piston_rad -pist_depth piston_rad -pist_depth 5 0 el_wl;
          piston_rad -pist_depth piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 5 0 el_wl
          -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 5 0 el_wl];

[xyb_low,topology_low]=nodegen(segments,'n');         % compute nodes and elements
M_low=size(xyb_low,1);N_low=size(topology_low,1);             % M nodes, N elementsq

% Excitation: Tweeter membrane displacement
nn_low=find(xyb_low(:,1)>=-piston_rad & xyb_low(:,1)<=piston_rad & xyb_low(:,2)>=-5e-4-eps-pist_depth & xyb_low(:,2)<= 5e-4-eps-pist_depth);

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];

% assign velocity vector
vn_low=zeros(M_low,1); vn_low(nn_low)=vampl* 1e-3;

theta_low=(0:fp_spc_deg*pi/180:thetamax)';rr_low=theta_low*Rmed;

fpxy_low=[Rmed*sin(theta_low) Rmed*cos(theta_low)];

Nfield_low = length(rr_low);
% -------------- High mesh density -------------

el_wl=6 * fr/c;
espac=1/el_wl;       
fp_spc_deg = asin(espac/Rmed)*180/pi;      
kp=2*pi*fr/c;        


segments=[-baffle_rad baffle_z -piston_rad baffle_z 5 0 el_wl;
-piston_rad baffle_z -piston_rad -pist_depth 5 0 el_wl;
          -piston_rad -pist_depth piston_rad -pist_depth 5 0 el_wl;
          piston_rad -pist_depth piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 5 0 el_wl
          -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 5 0 el_wl];

[xyb_high,topology_high]=nodegen(segments,'n');         % compute nodes and elements
M_high=size(xyb_high,1);N_high=size(topology_high,1);             % M nodes, N elementsq

% Excitation: Tweeter membrane displacement
nn_high=find(xyb_high(:,1)>=-piston_rad & xyb_high(:,1)<=piston_rad & xyb_high(:,2)>=-5e-4-eps-pist_depth & xyb_high(:,2)<= 5e-4-eps-pist_depth);

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];

% assign velocity vector
vn_high=zeros(M_high,1); vn_high(nn_high)=vampl* 1e-3;


theta_high=(0:fp_spc_deg*pi/180:thetamax)';rr_high=theta_high*Rmed;

fpxy_high=[Rmed*sin(theta_high) Rmed*cos(theta_high)];

Nfield_high = length(rr_high);
plot_BEM_mesh




%%



[A,B]=bem2d(xyb_low,topology_low,kp,betaP,xyb_chief);
    
disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])

% Pressure on the surface
B=1i*kp*rho*c*B;
ps=A\(-B*vn_low); % Solve the system
clear A B

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb_low,topology_low,kp,betaP,fpxy_low);
% solve the pressure on the field points
pF_low=(Ap*ps+1i*kp*rho*c*Bp*vn_low)./CConst;



[A,B]=bem2d(xyb_high,topology_high,kp,betaP,xyb_chief);
    
disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])

% Pressure on the surface
B=1i*kp*rho*c*B;
ps=A\(-B*vn_high); % Solve the system
clear A B

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb_high,topology_high,kp,betaP,fpxy_high);
% solve the pressure on the field points
pF_high=(Ap*ps+1i*kp*rho*c*Bp*vn_high)./CConst;

% dir_comp_mesh_size



% -------------- Low mesh density -------------
fr = 8000;
el_wl=1 * fr/c;
espac=1/el_wl;       
fp_spc_deg = asin(espac/Rmed)*180/pi;      
kp=2*pi*fr/c;        


segments=[-baffle_rad baffle_z -piston_rad baffle_z 5 0 el_wl;
-piston_rad baffle_z -piston_rad -pist_depth 5 0 el_wl;
          -piston_rad -pist_depth piston_rad -pist_depth 5 0 el_wl;
          piston_rad -pist_depth piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 5 0 el_wl
          -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 5 0 el_wl];

[xyb_low,topology_low]=nodegen(segments,'n');         % compute nodes and elements
M_low=size(xyb_low,1);N_low=size(topology_low,1);             % M nodes, N elementsq

% Excitation: Tweeter membrane displacement
nn_low=find(xyb_low(:,1)>=-piston_rad & xyb_low(:,1)<=piston_rad & xyb_low(:,2)>=-5e-4-eps-pist_depth & xyb_low(:,2)<= 5e-4-eps-pist_depth);

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];

% assign velocity vector
vn_low=zeros(M_low,1); vn_low(nn_low)=vampl* 1e-3;

theta_low2=(0:fp_spc_deg*pi/180:thetamax)';rr_low2=theta_low2*Rmed;

fpxy_low=[Rmed*sin(theta_low2) Rmed*cos(theta_low2)];

Nfield_low = length(rr_low2);
% -------------- High mesh density -------------

el_wl=6 * fr/c;
espac=1/el_wl;       
fp_spc_deg = asin(espac/Rmed)*180/pi;      
kp=2*pi*fr/c;        


segments=[-baffle_rad baffle_z -piston_rad baffle_z 5 0 el_wl;
-piston_rad baffle_z -piston_rad -pist_depth 5 0 el_wl;
          -piston_rad -pist_depth piston_rad -pist_depth 5 0 el_wl;
          piston_rad -pist_depth piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 5 0 el_wl
          -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 5 0 el_wl];

[xyb_high,topology_high]=nodegen(segments,'n');         % compute nodes and elements
M_high=size(xyb_high,1);N_high=size(topology_high,1);             % M nodes, N elementsq

% Excitation: Tweeter membrane displacement
nn_high=find(xyb_high(:,1)>=-piston_rad & xyb_high(:,1)<=piston_rad & xyb_high(:,2)>=-5e-4-eps-pist_depth & xyb_high(:,2)<= 5e-4-eps-pist_depth);

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];

% assign velocity vector
vn_high=zeros(M_high,1); vn_high(nn_high)=vampl* 1e-3;


theta_high2=(0:fp_spc_deg*pi/180:thetamax)';rr_high2=theta_high2*Rmed;

fpxy_high=[Rmed*sin(theta_high2) Rmed*cos(theta_high2)];

Nfield_high = length(rr_high2);




%%



[A,B]=bem2d(xyb_low,topology_low,kp,betaP,xyb_chief);
    
disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])

% Pressure on the surface
B=1i*kp*rho*c*B;
ps=A\(-B*vn_low); % Solve the system
clear A B

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb_low,topology_low,kp,betaP,fpxy_low);
% solve the pressure on the field points
pF_low2=(Ap*ps+1i*kp*rho*c*Bp*vn_low)./CConst;



[A,B]=bem2d(xyb_high,topology_high,kp,betaP,xyb_chief);
    
disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])

% Pressure on the surface
B=1i*kp*rho*c*B;
ps=A\(-B*vn_high); % Solve the system
clear A B

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb_high,topology_high,kp,betaP,fpxy_high);
% solve the pressure on the field points
pF_high2=(Ap*ps+1i*kp*rho*c*Bp*vn_high)./CConst;

dir_comp_mesh_size
