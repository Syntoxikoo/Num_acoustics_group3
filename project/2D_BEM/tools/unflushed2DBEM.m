function [fr,p_fieldUF,rr,theta] = unflushed2DBEM(el_wl,fr)
% 2D BEM unflushed
addpath(genpath("2DBEM"))


% -------------------- INPUT DATA ---------------------
 
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
espac=1/el_wl;       % Spacing of field points
fp_spc_deg = 1;      % Field point spacing in degree for the arc
kp=2*pi*fr/c;        % Wavenumber (1/m)
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

segments=[-baffle_rad baffle_z -piston_rad baffle_z 5 0 el_wl;
          -piston_rad baffle_z -piston_rad -pist_depth 5 0 el_wl;
          -piston_rad -pist_depth piston_rad -pist_depth 5 0 el_wl;
          piston_rad -pist_depth piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) -baffle_rad (baffle_z-thickness) 5 0 el_wl
          -baffle_rad (baffle_z-thickness) -baffle_rad baffle_z 5 0 el_wl];

[xyb,topology]=nodegen(segments,'n');         % compute nodes and elements
M=size(xyb,1);N=size(topology,1);             % M nodes, N elements

% CHIEF points:
xyb_chief_left=[linspace(0,-baffle_rad+0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief_right=[linspace(0,baffle_rad-0.05,5)', linspace(-pist_depth-0.05,baffle_z-thickness+0.05,5)', -ones(5,1)];
xyb_chief= [xyb_chief_left; xyb_chief_right];
%hold on; plot(xyb_chief(:,1),xyb_chief(:,2),'md',DisplayName="Chief")


% Excitation: Tweeter membrane displacement
nn1=find(xyb(:,1)>=-piston_rad & xyb(:,1)<=piston_rad & xyb(:,2)>=-5e-3-eps-pist_depth & xyb(:,2)<= 5e-3-eps-pist_depth);

%plot(xyb(nn1,1),xyb(nn1,2), "^b")

% assign velocity vector
vn=zeros(M,1); vn(nn1)=vampl* 1e-3;




% Definition of field points
% field points on the arc
theta=(0:fp_spc_deg*pi/180:thetamax)';rr=theta*Rmed;
fpxyC=[Rmed*sin(theta) Rmed*cos(theta)];
% field points on the radial direction
rfp=(Rmed*1.1:-espac:Rmed*0.9)';
fpxyR=[zeros(length(rfp),1) rfp];
fpxy=[fpxyC ; fpxyR];

% test figures showing velocity: z-component and normal velocity, and field points

%plot(fpxy(:,1),fpxy(:,2),'g+', DisplayName="Field points");
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);

%normplot=mean(mean(abs(diff(xyb(:,1:2))))); % quiver on the previous geometry figure
% quiver(rzb(:,1),rzb(:,2),nvect(:,1).*vn*normplot,nvect(:,2)*normplot.*vn,0,'-r');

%%

p_fieldUF = zeros(length(fpxy),length(fr));
for ii= 1:length(fr)
% CALCULATION OF RESULTS
    % BEM matrices calculation
    disp(['Calculating f= ' num2str(fr(ii)) ' Hz'])
    [A,B]=bem2d(xyb,topology,kp(ii),betaP,xyb_chief);
    
    disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])
    
    % Pressure on the surface
    B=1i*kp(ii)*rho*c*B;
    ps=A\(-B*vn); % Solve the system
    clear A B
    
    
    % figure;
    % subplot(2,1,1);
    % title("pressure at the surface, f ="+fr+"Hz")
    % plot((nn1),abs(ps(nn1)), "DisplayName","Ps_piston" );grid; hold on;
    % plot((max(nn1):length(ps)),abs(ps(max(nn1):end)), "DisplayName","Ps_baffle" )
    % xlabel('Nodes on the generator'); ylabel('Pressure modulus (Pa)')
    % subplot(2,1,2);plot(angle(ps)*180/pi);grid
    % xlabel('Nodes on the generator'); ylabel('Phase of the pressure (degrees)')
    % legend()

    % calculate corresponding rows of coefficients
    [Ap,Bp,CConst]=fieldpoints(xyb,topology,kp(ii),betaP,fpxy);
    % solve the pressure on the field points
    p_fieldUF(:,ii)=(Ap*ps+1i*kp(ii)*rho*c*Bp*vn)./CConst;
end
end