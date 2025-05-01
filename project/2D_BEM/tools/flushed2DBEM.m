function [fr,p_fieldF,rr,fp_theta] = flushed2DBEM(el_wl,fr,fp_theta)
% el_wl = element mesh density
% fr = frequency
% fp_theta = field point position in degree

% 2D BEM flushed
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
%fr=500;
%fr = [1000 2000 4000 8000];
vampl=1;             % Amplitude of the diaphragm movement (m/s) 
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
%piston_ctr = 0;
baffle_rad = 0.3;
baffle_z = 0;
thickness = 0.3; % changing thickness seems to have no effect since interior pb is omitted

segments=[-baffle_rad baffle_z -piston_rad baffle_z 1 0 el_wl;
          -piston_rad baffle_z piston_rad baffle_z 1 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 1 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 1 0 el_wl;
          baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 1 0 el_wl
          -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 1 0 el_wl];

[xyb,topology]=nodegen(segments,'n');         % compute nodes and elements
M=size(xyb,1); N=size(topology,1);             % M nodes, N elements

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(baffle_z-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(baffle_z-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];
%hold on; plot(xyb_chief(:,1),xyb_chief(:,2),'md',DisplayName="Chief")


% Excitation: Tweeter membrane displacement
nn1=find(xyb(:,1)>=-piston_rad-1e-5 & xyb(:,1)<=piston_rad+1e-5 & xyb(:,2)>=-1e-5 & xyb(:,2)<= 1e-5);

% plot(xyb(nn1,1),xyb(nn1,2), "^b")

% assign velocity vector
vn=zeros(M,1); vn(nn1)=vampl* 1e-3;




% Definition of field points
% field points on the arc
%fp_theta=(0:fp_spc_deg*pi/180:thetamax)';
rr=fp_theta*Rmed;
fpxyC=[Rmed*sin(fp_theta) Rmed*cos(fp_theta)];
% field points on the radial direction
% rfp=(Rmed*1.1:-espac:Rmed*0.9)';
% fpxyR=[zeros(length(rfp),1) rfp];
fpxy= fpxyC;%[fpxyC ; fpxyR];


%%

p_fieldF = zeros(length(fpxy),length(fr));
for ii= 1:length(fr)
% CALCULATION OF RESULTS
    % BEM matrices calculation
    %disp(['Calculating f= ' num2str(fr(ii)) ' Hz'])
    [A,B]=bem2d(xyb,topology,kp(ii),betaP,xyb_chief);
    
    %disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])
    
    % Pressure on the surface
    B=1i*kp(ii)*rho*c*B;
    ps=A\(-B*vn); % Solve the system
    clear A B

    % calculate corresponding rows of coefficients
    [Ap,Bp,CConst]=fieldpoints(xyb,topology,kp(ii),betaP,fpxy);
    % solve the pressure on the field points
    p_fieldF(:,ii)=(Ap*ps+1i*kp(ii)*rho*c*Bp*vn)./CConst;
end

end