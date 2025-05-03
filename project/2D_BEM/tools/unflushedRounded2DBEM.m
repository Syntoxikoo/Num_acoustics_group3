function [fr,p_field,rr,fp_theta,xyb,topology] = unflushedRounded2DBEM(el_wl,fr,fp_theta)
% el_wl = element mesh density
% fr = frequency
% fp_theta = field point position in degree

% 2D BEM unflushed
addpath(genpath("2DBEM"))


% -------------------- INPUT DATA ---------------------
 
% calculate radiant
fp_theta = fp_theta*pi/180;


% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 
% General parameters:
% fr=500;
% fr = [1000 2000 4000 8000];
vampl=1;             % Amplitude of the diaphragm movement (m/s) 
%el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
% espac=1/el_wl;       % Spacing of field points
% fp_spc_deg = 1;      % Field point spacing in degree for the arc
kp=2*pi*fr/c;        % Wavenumber (1/m)
betaP = NaN;           % normalised admittance of the plane, at k.


% Field points parameters
%thetamax=2*pi;        % Maximum angle for representation
Rmed=1;               % Radius of the arc of field points
%beta=0; 



% DEFINITION OF DOMAIN GEOMETRY AND EXCITATION
piston_rad = 0.1;
pist_depth = 0.1;
piston_ctr = 0;
baffle_rad = 0.3;
baffle_z = 0;
baffle_y = 0;
thickness = 0.3; % changing thickness seems to have no effect since interior pb is omitted

    segments=[-baffle_rad baffle_y -2*piston_rad baffle_y 1 0 el_wl;
        -2*piston_rad baffle_y -piston_rad -pist_depth 1 0.1 el_wl;
        -piston_rad -pist_depth piston_rad -pist_depth 1 0 el_wl;
        piston_rad -pist_depth 2*piston_rad baffle_y 5 0.1 el_wl;
        2*piston_rad baffle_y baffle_rad baffle_y 1 0 el_wl;
        baffle_rad baffle_y baffle_rad (baffle_y - thickness) 1 0 el_wl;
        baffle_rad (baffle_y-thickness) -baffle_rad (baffle_y-thickness) 1 0 el_wl
        -baffle_rad (baffle_y-thickness) -baffle_rad baffle_y 1 0 el_wl];

[xyb,topology]=nodegen(segments,'y');         % compute nodes and elements
M=size(xyb,1);N=size(topology,1);             % M nodes, N elements

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];
hold on; plot(xyb_chief(:,1),xyb_chief(:,2),'md',DisplayName="Chief")


% Excitation: Tweeter membrane displacement
nn1=find(xyb(:,1)>=-piston_rad-1e-5 & xyb(:,1)<=piston_rad + 1e-5 & xyb(:,2)>= -1e-5-pist_depth & xyb(:,2)<= 1e-5-pist_depth);

plot(xyb(nn1,1),xyb(nn1,2), "^b")

% assign velocity vector
vn=zeros(M,1); vn(nn1)=vampl* 1e-3;




% Definition of field points
% field points on the arc
% fp_theta=(0:fp_spc_deg*pi/180:thetamax)';
rr=fp_theta*Rmed;
fpxyC=[Rmed*sin(fp_theta) Rmed*cos(fp_theta)];
% field points on the radial direction
%rfp=(Rmed*1.1:-espac:Rmed*0.9)';
%fpxyR=[zeros(length(rfp),1) rfp];
fpxy=fpxyC;%[fpxyC ; fpxyR];

% test figures showing velocity: z-component and normal velocity, and field points

plot(fpxy(:,1),fpxy(:,2),'b*', DisplayName="Field points");
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

% normplot=mean(mean(abs(diff(xyb(:,1:2))))); % quiver on the previous geometry figure
% quiver(rzb(:,1),rzb(:,2),nvect(:,1).*vn*normplot,nvect(:,2)*normplot.*vn,0,'-r');

%%

p_field = zeros(length(fpxy),length(fr));
for ii= 1:length(fr)
% CALCULATION OF RESULTS
    % BEM matrices calculation
    disp(['Calculating f= ' num2str(fr(ii)) ' Hz'])
    [A,B]=bem2d(xyb,topology,kp(ii),betaP,xyb_chief);
    
    %disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])
    
    % Pressure on the surface
    B=1i*kp(ii)*rho*c*B;
    ps=A\(-B*vn); % Solve the system
    clear A B

    % calculate corresponding rows of coefficients
    [Ap,Bp,CConst]=fieldpoints(xyb,topology,kp(ii),betaP,fpxy);
    % solve the pressure on the field points
    p_field(:,ii)=(Ap*ps+1i*kp(ii)*rho*c*Bp*vn)./CConst;
end
end