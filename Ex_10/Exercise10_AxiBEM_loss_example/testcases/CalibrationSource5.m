
% Sound source for secondary calibration of microphones by comparison.
% Joint project with DFM.

% Version 5: Uses a Vifa tweeter in a cut "tennis ball". This is the
% existing source used for calibration by substitution.

% The acoustic center is calculated from field points.

% Vicente Cutanda Henríquez 05-2007
% Johan Heilmann 04-2007


clear
% close all
 

% INPUT DATA
 
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
 
% General parameters:
fr=1000;             % Frequency (Hz)    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<**************
vampl=1;             % Amplitude of the diaphragm movement (m/s) 
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points
kp=2*pi*fr/c;        % Wavenumber (1/m)
m=0;                 % Axisymetrical excitation, circunferential mode m=0

% Source's geometrical parameters (m):
Rt=34e-3;            % Radius of the tweeter
Rcmax=0.015;         % Largest radius of the cone base
Rcmin=0.008;         % Shortest radius of the cone base
Hc=2e-2;             % Heigth of the cone base
Hr=5e-3;             % Thickness of the Tweeter rim

% Field points parameters
thetamax=pi/6;        % Maximum angle for representation
Rmed=1;               % Radius of the arc of field points
beta=0; % beta=pi/18; % Angle at wich calculate the radial variation



% DEFINITION OF DOMAIN GEOMETRY AND EXCITATION

[tweetersegments,m_min,m_mid,m_max,tweeter_end] = tweeter(c/max(fr),0);

segments=[tweetersegments;
          Rt 0 Rt -Hr 1 Rt el_wl;
          Rt -Hr Rcmax -sqrt(Rt^2-Rcmax^2)-Hr 10 Rt el_wl;
          Rcmax -sqrt(Rt^2-Rcmax^2)-Hr Rcmin -Rt-Hc-Hr 6 0 el_wl;
          Rcmin -Rt-Hc-Hr 0 -Rcmin-Rt-Hc-Hr 5 Rcmin el_wl];

[rzb,topology]=nodegen(segments,'n');         % compute nodes and elements
[nvect]=normals(rzb,topology,'y');            % normal vectors
M=size(rzb,1);N=size(topology,1);             % M nodes, N elements

% CHIEF points:
rzb_chief=[linspace(0,0.02,10)' linspace(-0.02,-0.01,10)' -ones(10,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md',DisplayName="Chief")

% Excitation: Tweeter membrane displacement
vz=zeros(M,1);
nn1=find(rzb(:,1)>=m_min(1) & rzb(:,1)<=m_mid(1) & rzb(:,2)>=-5e-3-eps);
nn2=find(rzb(:,1)>=m_mid(1) & rzb(:,1)<=m_max(1) & rzb(:,2)>=-5e-3-eps);

vz(nn1,1)=(rzb(nn1,1)-m_min(1))*vampl/(m_mid(1)-m_min(1));
vz(nn2,1)=-(rzb(nn2,1)-m_mid(1))*vampl/(m_max(1)-m_mid(1)) + vampl;

vn=dot([zeros(M,1) vz]',nvect')';

% Definition of field points
% field points on the arc
theta=(0:2*asin(espac/2/Rmed):thetamax)';rr=theta*Rmed;
fprzC=[Rmed*sin(theta) Rmed*cos(theta) ones(length(theta),1)*4*pi];
% field points on the radial direction
rfp=(Rmed*1.1:-espac:Rmed*0.9)';
fprzR=[rfp*sin(beta) rfp*cos(beta) ones(length(rfp),1)*4*pi];
fprz=[fprzC ; fprzR];

% test figures showing velocity: z-component and normal velocity, and field points
hold on 
leg1 = plot(m_mid(1),m_mid(2),"b^","markersize",10, "DisplayName","mid");
leg2 = plot(m_max(1),m_max(2),"r^","markersize",10,"DisplayName","max");
leg3 = plot(m_min(1),m_min(2),"m^","markersize",10,"DisplayName","min");
plot(fprz(:,1),fprz(:,2),'g+', DisplayName="Field points");
legend([leg1,leg2,leg3])
normplot=mean(mean(abs(diff(rzb(:,1:2))))); % quiver on the previous geometry figure
quiver(rzb(:,1),rzb(:,2),nvect(:,1).*vn*normplot,nvect(:,2)*normplot.*vn,0,'-r');


figure;plot(rzb([nn1;nn2],1),vn([nn1;nn2]),"DisplayName", "normal vel");hold on
    plot(rzb([nn1;nn2],1),vz([nn1;nn2]),"DisplayName","vel in z dir");grid; 
plot(rzb([1;nn1],1),vz([1;nn1]),"DisplayName","vel z nn0")% velocity shape figure
% plot(rzb([nn2;end-1],1),vz([nn2;end-1]),"DisplayName","vel z nn3")
title("velocity shape fig")
xlabel("rho")
ylabel("vel")
legend

% CALCULATION OF RESULTS

% BEM matrices calculation
disp(['Calculating f= ' num2str(fr) ' Hz'])
[A,B,CConst]=BEMEquat0(rzb,topology,kp,m,rzb_chief);
disp(['Condition numbers, A: ' num2str(cond(A)) ' , B: ' num2str(cond(B))])

% Pressure on the surface
B=i*kp*rho*c*B; % B times i*k*rho*c to compute pressure instead of velocity potential
ps=A\(-B*vn); % Solve the system
clear A B

% Figure to test the surface solution
figure;
subplot(2,1,1);plot(abs(ps));grid
xlabel('Nodes on the generator'); ylabel('Pressure modulus (Pa)')
subplot(2,1,2);plot(angle(ps)*180/pi);grid
xlabel('Nodes on the generator'); ylabel('Phase of the pressure (degrees)')

% Field points calculation
p_field=FieldPnt3(fprz,ps,vn,rzb,topology,kp,m,rho,c);

% Acoustc centre estimation
p_far=abs(p_field(length(rr)+1));
p_near=abs(p_field(length(rr)+length(rfp)));
r_AcCen=rfp(1) - (rfp(1)-rfp(end))/(1/p_far-1/p_near)/p_far;


% PRESENTATION OF RESULTS

% Sound pressure on an arc of far field points
figure; 
subplot(3,2,1);plot(rr,20*log10(abs(p_field(1:length(rr)))/20e-6));grid
xlabel(['Arc length, for R = ' num2str(Rmed) ' m']); ylabel('SPL [dB]')
title(['p modulus (dB), freq.= ' num2str(fr) ' Hz']);
subplot(3,2,3);plot(rr,angle(p_field(1:length(rr)))*180/pi);grid
xlabel(['Arc length, for R = ' num2str(Rmed) ' m']); ylabel('Phase of the pressure [degrees]')
subplot(3,2,5);plot(rr,gradient(20*log10(abs(p_field(1:length(rr)))),espac)/100,'-x');grid
xlabel(['Arc length, for R = ' num2str(Rmed) ' m']); ylabel('Gradient of the pressure [dB/cm]')


% Sound pressure on an radius of far field points
subplot(3,2,2);plot(rfp,20*log10(abs(p_field(length(rr)+1:length(rr)+length(rfp)))/20e-6));grid
xlabel(['Distance on the axis z (m)']); ylabel('SPL [dB]')
title(['p modulus (dB), freq.= ' num2str(fr) ' Hz']);
subplot(3,2,4);plot(rfp,angle(p_field(length(rr)+1:length(rr)+length(rfp)))*180/pi);grid
xlabel(['Distance on the axis z (m)']); ylabel('Phase of the pressure [degrees]')
subplot(3,2,6);plot(rfp,gradient(20*log10(abs(p_field(length(rr)+1:length(rr)+length(rfp)))),espac)/100,'-x');grid
xlabel(['Distance on the axis z (m)']); ylabel('Gradient of the pressure [dB/cm]')

