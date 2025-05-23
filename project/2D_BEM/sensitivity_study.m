clear; close all; 

% Problem Definition

% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c0,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 
% General parameters:

fr =  linspace(1500,4000,100);
vampl=1;             % Amplitude of the diaphragm movement (m/s) 
betaP = NaN;           % normalised admittance of the plane, at k.
% Field points parameters
theta = deg2rad(60);
beta=0; 
Rmed = 1;

% DEFINITION OF DOMAIN GEOMETRY AND EXCITATION
piston_rad = 0.1;
pist_depth = 0.1;
piston_ctr = 0;
baffle_rad = 0.3;
baffle_z = 0;
thickness = 0.3; % changing thickness seems to have no effect since interior pb is omitted

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];

rr=theta*Rmed;

fpxy=[Rmed*sin(theta) Rmed*cos(theta)];

pF_1 = zeros(length(fr),1);
for ii = 1: length(fr)
    disp("Solving problem unflushed : "+ii+"/"+length(fr))
    el_wl=6 * (fr(ii))/c0;
    espac=1/el_wl;       
        
    kp=2*pi*(fr(ii))/c0;        


    segments=[-baffle_rad baffle_z -piston_rad baffle_z 5 0 el_wl;
    -piston_rad baffle_z -piston_rad -pist_depth 5 0 el_wl;
            -piston_rad -pist_depth piston_rad -pist_depth 5 0 el_wl;
            piston_rad -pist_depth piston_rad baffle_z 5 0 el_wl;
            piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
            baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
            baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 5 0 el_wl
            -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 5 0 el_wl];

    [xyb,topology]=nodegen(segments,'n');         % compute nodes and elements
    M=size(xyb,1);N=size(topology,1);             % M nodes, N elementsq

    % Excitation: Tweeter membrane displacement
    nn=find(xyb(:,1)>=-piston_rad & xyb(:,1)<=piston_rad & xyb(:,2)>=-5e-4-eps-pist_depth & xyb(:,2)<= 5e-4-eps-pist_depth);

    % assign velocity vector
    vn=zeros(M,1); vn(nn)=vampl* 1e-3;


    % --------- Solving -----------
    [A,B]=bem2d(xyb,topology,kp,betaP,xyb_chief);
    disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])

    % Pressure on the surface
    B=1i*kp*rho*c0*B;
    ps=A\(-B*vn); % Solve the system
    clear A B
    % calculate corresponding rows of coefficients
    [Ap,Bp,CConst]=fieldpoints(xyb,topology,kp,betaP,fpxy);
    % solve the pressure on the field points
    pF_1(ii)=(Ap*ps+1i*kp*rho*c0*Bp*vn)./CConst;
end


pF_2 = zeros(length(fr),1);
for ii = 1: length(fr)
    disp("Solving problem horn : "+ii+"/"+length(fr))
    el_wl=6 * (fr(ii))/c0;
    espac=1/el_wl;       
        
    kp = 2*pi*(fr(ii))/c0;        


    % DEFINITION OF DOMAIN GEOMETRY AND EXCITATION
    piston_rad = 0.1;
    pist_depth = 0.1;
    piston_ctr = 0;
    baffle_rad = 0.3;
    baffle_y = 0;
    thickness = 0.3; % changing thickness seems to have no effect since interior pb is omitted

    segments=[-baffle_rad baffle_y -2*piston_rad baffle_y 1 0 el_wl;
    -2*piston_rad baffle_y -piston_rad -pist_depth 1 0.1 el_wl;
    -piston_rad -pist_depth piston_rad -pist_depth 1 0 el_wl;
    piston_rad -pist_depth 2*piston_rad baffle_y 1 0.1 el_wl;
    2*piston_rad baffle_y baffle_rad baffle_y 1 0 el_wl;
    baffle_rad baffle_y baffle_rad (baffle_y - thickness) 1 0 el_wl;
    baffle_rad (baffle_y-thickness) -baffle_rad (baffle_y-thickness) 1 0 el_wl
    -baffle_rad (baffle_y-thickness) -baffle_rad baffle_y 1 0 el_wl];
    
    [xyb2,topology2]=nodegen(segments,'n');         % compute nodes and elements
    M=size(xyb2,1);N=size(topology2,1);             % M nodes, N elementsq

    % Excitation: Tweeter membrane displacement
    nn2=find(xyb2(:,1)>=-piston_rad-1e-5 & xyb2(:,1)<=piston_rad + 1e-5 & xyb2(:,2)>= -1e-5-pist_depth & xyb2(:,2)<= 1e-5-pist_depth);

    % assign velocity vector
    vn=zeros(M,1); vn(nn2)=vampl* 1e-3;


    % --------- Solving -----------
    [A,B]=bem2d(xyb2,topology2,kp,betaP,xyb_chief);
    disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])

    % Pressure on the surface
    B=1i*kp*rho*c0*B;
    ps=A\(-B*vn); % Solve the system
    clear A B
    % calculate corresponding rows of coefficients
    [Ap,Bp,CConst]=fieldpoints(xyb2,topology2,kp,betaP,fpxy);
    % solve the pressure on the field points
    pF_2(ii)=(Ap*ps+1i*kp*rho*c0*Bp*vn)./CConst;
end
save("sensitivity_study_recessed_vs_bevelcurve2","fr","pF_1","pF_2","theta","fpxy")
%% 
plot_sensitivity_study_Bem


