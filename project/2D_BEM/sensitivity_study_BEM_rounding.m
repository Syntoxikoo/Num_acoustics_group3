addpath(genpath("tools"))
addpath(genpath("2DBEM"))
clear;
clc;
close all;

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



rr=theta*Rmed;

fpxy=[Rmed*sin(theta) Rmed*cos(theta)];
%% first one
pF_1 = zeros(length(fr),1);
for ii = 1: length(fr)
    disp("Solving problem unflushed : "+ii+"/"+length(fr))
    el_wl=6 * (fr(ii))/c0;
    espac=1/el_wl;       
        
    kp=2*pi*(fr(ii))/c0;        


    segments=[-baffle_rad baffle_z -piston_rad baffle_z 1 0 el_wl;
    -piston_rad baffle_z -piston_rad -pist_depth 1 0 el_wl;
            -piston_rad -pist_depth piston_rad -pist_depth 1 0 el_wl;
            piston_rad -pist_depth piston_rad baffle_z 1 0 el_wl;
            piston_rad baffle_z baffle_rad baffle_z 1 0 el_wl;
            baffle_rad baffle_z baffle_rad (baffle_z - thickness) 1 0 el_wl;
            baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 1 0 el_wl
            -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 1 0 el_wl];

    [xyb,topology]=nodegen(segments,'n');         % compute nodes and elements
    M=size(xyb,1);N=size(topology,1);             % M nodes, N elementsq

    % CHIEF points:
    xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
    xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
    xyb_chief= [xyb_chief_left; xyb_chief_right];
    % hold on; plot(xyb_chief(:,1),xyb_chief(:,2),'md',DisplayName="Chief")
    % Excitation: Tweeter membrane displacement
    nn1=find(xyb(:,1)>=-piston_rad-1e-5 & xyb(:,1)<=piston_rad + 1e-5 & xyb(:,2)>= -1e-5-pist_depth & xyb(:,2)<= 1e-5-pist_depth);
    % plot(xyb(nn1,1),xyb(nn1,2), "^b")
    % assign velocity vector
    vn=zeros(M,1); vn(nn1)=vampl* 1e-3;


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

%% second one
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

    % CHIEF points:
    xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
    xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
    xyb_chief= [xyb_chief_left; xyb_chief_right];
    % hold on; plot(xyb_chief(:,1),xyb_chief(:,2),'md',DisplayName="Chief")
    % Excitation: Tweeter membrane displacement
    nn2=find(xyb2(:,1)>=-piston_rad-1e-5 & xyb2(:,1)<=piston_rad + 1e-5 & xyb2(:,2)>= -1e-5-pist_depth & xyb2(:,2)<= 1e-5-pist_depth);
    % plot(xyb2(nn2,1),xyb2(nn2,2), "^b")
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
save("data/sensitivity_study_rounded","fr","pF_1","pF_2","theta","fpxy","xyb","xyb2","topology","topology2")
%% Plot amplitude
nrows = 1;
ncols = 1;
heightScale = 0.7; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

tile1 = nexttile;

Leg(1) = semilogx(fr, 20 * log10(abs(pF_1)/2e-5), 'LineWidth', 1.5,Color="#181748"); hold on ; grid on 
Leg(2) = semilogx(fr, 20 * log10(abs(pF_2)/2e-5), 'LineWidth', 1.5,Color="#810100");

xlabel("Frequency (Hz)")
ylabel("Level (dB re 20 \mu Pa)")
ylim([40 80])
leg1 = legend(Leg, "Recessed","Recessed and rounded",NumColumns=2);
leg1.Layout.Tile = 'north';
title("Effect of change in piston depth " +rad2deg(theta)+"° from the axis")
hold off
saveas(gcf,'figures/sensitivity_study_BEM_rounded.svg')

%% Plot Setup
nrows = 1;
ncols = 2;
heightScale = 0.55; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

tile1 = nexttile;

for ii =1:size(topology,1)
    Leg(1) = plot(xyb(topology(ii,1:3),1),xyb(topology(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on 
end
Leg(2) = plot(fpxy(1),fpxy(2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);

xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

grid on ; hold off
title( "Recessed 10 cm")


tile2 = nexttile;

for ii =1:size(topology2,1)
    plot(xyb2(topology2(ii,1:3),1),xyb2(topology2(ii,1:3),2),"Color","black","Linewidth",1.5); 
    hold on;
end
plot(fpxy(1),fpxy(2),"Color", "#810100","Marker",".","LineStyle","none", "MarkerSize",10);
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);



xlabel(tiled, "x (m)")
ylabel(tiled, "y (m)")
title( "Recessed and rounded cm")
grid on;
hold off;

leg = legend(Leg, {"Meshed body", "Observation point"}, 'NumColumns', 3); 
leg.Layout.Tile = 'north'; 

saveas(gcf,'figures/geometry_sensitivity_study_rounded.svg')



