clear


% -------------------- INPUT DATA ---------------------
 
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 
% General parameters:
fr=500;
fr = [1000 2000 4000 8000];
vampl=1;             % Amplitude of the diaphragm movement (m/s) 
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points
fp_spc_deg = 1;      % Field point spacing in degree for the arc
kp=2*pi*fr/c;        % Wavenumber (1/m)
m=0;                 % Axisymetrical excitation, circunferential mode m=0


% Field points parameters
thetamax=pi;        % Maximum angle for representation
Rmed=1;               % Radius of the arc of field points
beta=0; 



% DEFINITION OF DOMAIN GEOMETRY AND EXCITATION
pist_depth = 0.1;
piston_rad = 0.1;
piston_ctr = 0;
baffle_rad = 0.3;
baffle_z = 0;
thickness = 0.3; % changing thickness seems to have no effect since interior pb is omitted

segments=[0 -pist_depth piston_rad -pist_depth 5 0 el_wl;
          piston_rad -pist_depth piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) 0 (baffle_z-thickness) 5 0 el_wl];

[rzb,topology]=nodegen(segments,'n');         % compute nodes and elements
[nvect]=normals(rzb,topology,'y');            % normal vectors
M=size(rzb,1);N=size(topology,1);             % M nodes, N elements

% CHIEF points:
rzb_chief=[linspace(0,baffle_rad-0.05,10)' linspace(baffle_z-pist_depth-0.05,baffle_z-thickness+0.05,10)' -ones(10,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md',DisplayName="Chief")


% Excitation: Tweeter membrane displacement
vz=zeros(M,1);
nn1=find(rzb(:,1)>=piston_ctr & rzb(:,1)<=piston_rad & rzb(:,2)>=-5e-3-eps-pist_depth & rzb(:,2)<= 5e-3-eps-pist_depth);

plot(rzb(nn1,1),rzb(nn1,2), "^b")

% vz(nn1,1)=(rzb(nn1,1)-m_min(1))*vampl/(m_mid(1)-m_min(1));
vz(nn1,1) = vampl* 1e-3;

vn=dot([zeros(M,1) vz]',nvect')';



% Definition of field points
% field points on the arc
theta=(0:fp_spc_deg*pi/180:thetamax)';rr=theta*Rmed;
fprzC=[Rmed*sin(theta) Rmed*cos(theta) ones(length(theta),1)*4*pi];
% field points on the radial direction
rfp=(Rmed*1.1:-espac:Rmed*0.9)';
fprzR=[rfp*sin(beta) rfp*cos(beta) ones(length(rfp),1)*4*pi];
fprz=[fprzC ; fprzR];

% test figures showing velocity: z-component and normal velocity, and field points

plot(fprz(:,1),fprz(:,2),'g+', DisplayName="Field points");

normplot=mean(mean(abs(diff(rzb(:,1:2))))); % quiver on the previous geometry figure
% quiver(rzb(:,1),rzb(:,2),nvect(:,1).*vn*normplot,nvect(:,2)*normplot.*vn,0,'-r');

%%

p_fieldUF = zeros(length(fprz),length(fr));
for ii= 1:length(fr)
% CALCULATION OF RESULTS
    % BEM matrices calculation
    disp(['Calculating f= ' num2str(fr(ii)) ' Hz'])
    [A,B,CConst]=BEMEquat0(rzb,topology,kp(ii),m,rzb_chief,[],false);
    
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

    % Field points calculation
    p_fieldUF(:,ii)=FieldPnt3(fprz,ps,vn,rzb,topology,kp(ii),m,rho,c);
end

% Acoustc centre estimation
p_far=abs(p_fieldUF(length(rr)+1));
p_near=abs(p_fieldUF(length(rr)+length(rfp)));
r_AcCen=rfp(1) - (rfp(1)-rfp(end))/(1/p_far-1/p_near)/p_far;


% PRESENTATION OF RESULTS

% Sound pressure on an arc of far field points
% figure; 
% subplot(3,2,1);plot(rr,20*log10(abs(p_fieldUF(1:length(rr)))/20e-6));grid
% xlabel(['Arc length, for R = ' num2str(Rmed) ' m']); ylabel('SPL [dB]')
% title(['p modulus (dB), freq.= ' num2str(fr) ' Hz']);
% subplot(3,2,3);plot(rr,angle(p_fieldUF(1:length(rr)))*180/pi);grid
% xlabel(['Arc length, for R = ' num2str(Rmed) ' m']); ylabel('Phase of the pressure [degrees]')
% subplot(3,2,5);plot(rr,gradient(20*log10(abs(p_fieldUF(1:length(rr)))),espac)/100,'-x');grid
% xlabel(['Arc length, for R = ' num2str(Rmed) ' m']); ylabel('Gradient of the pressure [dB/cm]')
% 


% Sound pressure on an radius of far field points
% subplot(3,2,2);plot(rfp,20*log10(abs(p_fieldUF(length(rr)+1:length(rr)+length(rfp)))/20e-6));grid
% xlabel(['Distance on the axis z (m)']); ylabel('SPL [dB]')
% title(['p modulus (dB), freq.= ' num2str(fr) ' Hz']);
% subplot(3,2,4);plot(rfp,angle(p_fieldUF(length(rr)+1:length(rr)+length(rfp)))*180/pi);grid
% xlabel(['Distance on the axis z (m)']); ylabel('Phase of the pressure [degrees]')
% subplot(3,2,6);plot(rfp,gradient(20*log10(abs(p_fieldUF(length(rr)+1:length(rr)+length(rfp)))),espac)/100,'-x');grid
% xlabel(['Distance on the axis z (m)']); ylabel('Gradient of the pressure [dB/cm]')

%%
% figure;
% % Normalize the SPL values relative to the on-axis (0 degrees) response
% spl_values = 20*log10(abs(p_field(1:length(rr)))/20e-6);
% normalized_spl = spl_values - max(spl_values);
% 
% % Create polar plot with better formatting
% polarplot(theta, normalized_spl, 'LineWidth', 2);
% grid on;
% title(['Directivity Pattern at ' num2str(fr) ' Hz']);
% rlim([min(normalized_spl(:,1))-5, 0]);  % Set reasonable limits for radial axis
% rticks(-30:6:0);  % Set radial ticks every 6 dB
% thetaticks(0:15:90);  % Set angular ticks every 15 degrees
% thetalim([0 90]);  % Limit to the calculated range (0 to 90 degrees)


figure;

theta_full = [theta; pi + theta];


for i = 1:length(fr)
    p = p_fieldUF(:,i);
    spl_values = 20*log10(abs(p(1:length(rr)))/20e-6);
    normalized_spl = spl_values - max(spl_values);
    spl_full = [normalized_spl; flip(normalized_spl)];
    polarplot(theta_full, spl_full, 'LineWidth', 2, 'DisplayName', [num2str(fr(i)) ' Hz']);
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

save("data/BEM_uf","p_fieldUF","theta","fr","rr");