addpath(genpath("project"))
 
% -------------------- INPUT DATA ---------------------
 
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 
% General parameters:

fr = [1000 2000 4000 8000];
vampl=1;             % Amplitude of the diaphragm movement (m/s) 
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points
kp=2*pi*fr/c;        % Wavenumber (1/m)
m=0;                 % Axisymetrical excitation, circunferential mode m=0


% Field points parameters
thetamax=pi;        % Maximum angle for representation
Rmed=1;               % Radius of the arc of field points
fp_spc_deg = asin(espac/Rmed)*180/pi;
beta=0; 



% DEFINITION OF DOMAIN GEOMETRY AND EXCITATION
piston_rad = 0.1;
piston_ctr = 0;
baffle_rad = 0.3;
baffle_z = 0;
thickness = 0.3; % changing thickness seems to have no effect since interior pb is omitted



segments=[0 baffle_z piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) 0 (baffle_z-thickness) 5 0 el_wl];

[rzb,topology]=nodegen(segments,'n');         % compute nodes and elements
[nvect]=normals(rzb,topology,'y');            % normal vectors
M=size(rzb,1);N=size(topology,1);             % M nodes, N elements

% CHIEF points:
rzb_chief=[linspace(0,baffle_rad-0.05,10)' linspace(baffle_z-0.05,baffle_z-thickness+0.05,10)' -ones(10,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md',DisplayName="Chief")


% Excitation: Tweeter membrane displacement
vz=zeros(M,1);
nn1=find(rzb(:,1)>=piston_ctr & rzb(:,1)<=piston_rad & rzb(:,2)>=-5e-3-eps & rzb(:,2)<= 5e-3-eps);

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


p_fieldF_arr = cell(1,length(fr));
theta_arr = cell(1,length(fr));
rr_arr = cell(1, length(fr));
for ii= 1:length(fr)
    % Field definition
    el_wl=6 * fr(ii)/c;   
    espac = 1/el_wl;
    fp_spc_deg = asin(espac/Rmed)*180/pi;
    %Solve Geometry
    segments=[0 baffle_z piston_rad baffle_z 5 0 el_wl;
          piston_rad baffle_z baffle_rad baffle_z 5 0 el_wl;
          baffle_rad baffle_z baffle_rad (baffle_z - thickness) 5 0 el_wl;
          baffle_rad (baffle_z-thickness) 0 (baffle_z-thickness) 5 0 el_wl];

    [rzb,topology]=nodegen(segments,'n');         
    [nvect]=normals(rzb,topology,'n');            
    M=size(rzb,1);N=size(topology,1);             


    vz=zeros(M,1);
    nn1=find(rzb(:,1)>=piston_ctr & rzb(:,1)<=piston_rad & rzb(:,2)>=-5e-3-eps & rzb(:,2)<= 5e-3-eps);
    vz(nn1,1) = vampl* 1e-3;

    vn=dot([zeros(M,1) vz]',nvect')';

    theta=(0:fp_spc_deg*pi/180:thetamax)';rr=theta*Rmed;
    rr_arr{ii} = rr;
    theta_arr{ii} = theta;

    fprzC=[Rmed*sin(theta) Rmed*cos(theta) ones(length(theta),1)*4*pi];
    
    rfp=(Rmed*1.1:-espac:Rmed*0.9)';
    fprzR=[rfp*sin(beta) rfp*cos(beta) ones(length(rfp),1)*4*pi];
    fprz=[fprzC ; fprzR];

% CALCULATION OF RESULTS
    % BEM matrices calculation
    disp(['Calculating f= ' num2str(fr(ii)) ' Hz'])
    [A,B,CConst]=BEMEquat0(rzb,topology,kp(ii),m,rzb_chief,[],false);
    
    disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])
    
    % Pressure on the surface
    B=1i*kp(ii)*rho*c*B;
    ps=A\(-B*vn); % Solve the system
    clear A B
    
    
    % Field points calculation
    p_fieldF_arr{ii}=FieldPnt3(fprz,ps,vn,rzb,topology,kp(ii),m,rho,c);
end


%%
figure;




for ii = 1:length(fr)
    theta = cell2mat(theta_arr(ii));
    theta_full = [theta; pi + theta];
    p = cell2mat(p_fieldF_arr(ii));
    spl_values = 20*log10(abs(p(1:length(cell2mat(rr_arr(ii)))))/20e-6);
    normalized_spl = spl_values - max(spl_values);
    spl_full = [normalized_spl; flip(normalized_spl)];
    polarplot(theta_full, spl_full, 'LineWidth', 2, 'DisplayName', [num2str(fr(ii)) ' Hz']);
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

save("project/data/axi_BEM_f","p_fieldF_arr","theta_arr","fr","rr");