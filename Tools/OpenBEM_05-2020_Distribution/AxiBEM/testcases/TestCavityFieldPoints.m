% Example of sound field in a cavity with calculation of a mesh of field points
% A cylindrical cavity with a spherical cap
% 
% Vicente Cutanda 3-2008

clear
close all

% INPUT DATA

% Cavity dimensions (m):
Rc=0.05; % Cavity Radius
Hc=0.08; % Cavity Heigth
Hs=0.03; % Heigth og the spherical cap. If zero, the cavty is a cylinder
RTransd=0.0175;  % Radius of the vibrating piston (excitation)

% Other parameters:
fr=1000;          % Frequency (Hz)
espac=Rc/10+eps;  % Field point mesh spacing
vampl=1;          % Amplitude of the piston velocity (m/s) 
m=0;              % Axisymmetrical excitation, circumferential mode m=0

% physical constants (MKS units) of air in normal conditions:
pa = 101325;          % Static pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Relative humidity (%)
[rho,c]=amb2prop(pa,t,Hr,1000); % c: speed of sound (m/s); rho: air density (Kg/m3)
kp=2*pi*fr/c;     % Wavelength (1/m)


% DOMAIN GEOMETRY
% Minimum mesh density as a function of the highest frequency (see nodegen):
el_wl=6*max(fr)/c;   % Elements per meter, six elements per wavelength

% cavity generator
% The objects need to be defined in the clockwise direction by convention
if Hs==0 % cylinder case
    segments=[0 Hc Rc Hc 10 0 el_wl;...
              Rc Hc Rc 0 10 0 el_wl;Rc 0 0 0 10 0 el_wl];
else % cylinder with spherical cap
    RCurv=0.5*sqrt(Hs^2+Rc^2)/sin(atan(Hs/Rc)); % Curvature radius of the cap
    segments=[0 Hc+Hs Rc Hc 10 RCurv el_wl; ...
              Rc Hc Rc 0 10 0 el_wl;Rc 0 0 0 10 0 el_wl];
end


[rzb,topology]=nodegen(segments,'y'); % Nodes and elements computation
M=size(rzb,1);N=size(topology,1); % M nodes, N elements

% Vector of normal velocities
nn=find(rzb(:,1)<=RTransd & rzb(:,2)==0);
vs=zeros(M,1);vs(nn)=vampl;
hold on; plot(rzb(nn,1),rzb(nn,2),'r*')
% Set interior domain
rzb(:,end)=-rzb(:,end);
topology(:,end)=-topology(:,end);


% BEM MATRICES CALCULATION
disp(['Calculando f= ' num2str(fr) ' Hz'])
[A,B,CConst]=BEMEquat0(rzb,topology,kp,m);
disp(['condition numbers, of A: ' num2str(cond(A)) ' and of B: ' num2str(cond(B))])


% SOLUTION OVER THE SURFACE
% B is multiplied by i*k*rho*c to compute sound pressure instead of velocity potential
B=i*kp*rho*c*B;
%  solve system
ps=(A)\(-B*vs);

% Plot solution over the generator
figure;
subplot(2,1,1);plot(abs(ps));grid
xlabel('Nodes over the generator'); ylabel('Pressure modulus (Pa)')
subplot(2,1,2);plot(angle(ps)*180/pi);grid
xlabel('Nodes over the generator'); ylabel('Pressure phase (degrees)')


% FIELD POINT DEFINITION AND CALCULATION
rr=0:espac:max(rzb(:,1)); Nrr=length(rr);
zz=min(rzb(:,2)):espac:max(rzb(:,2)); Nzz=length(zz); 
[RR,ZZ]=meshgrid(rr,zz); % Mesh definition on the rho-z plane

%fprz=[RR(1:end)' ZZ(1:end)']; % in case FieldPnt3 is used
[p_field,RRb,ZZb]=FieldPnt(RR,ZZ,ps,vs,rzb,topology,kp,m,rho,c); % See FieldPnt help text


% PLOT SOUND PRESSURE ON THE FIELD POINT MESH
figure
subplot(1,2,1)
contourf(RRb,ZZb,abs(p_field))
title(['p modulus (Pa), freq.= ' num2str(fr) ' Hz']);xlabel('r (m)');ylabel('z (m)')
axis equal;axis([0 max(rzb(:,1)) min(rzb(:,2)) max(rzb(:,2))])
colorbar
subplot(1,2,2)
contourf(RRb,ZZb,angle(p_field)*180/pi)
title(['p phase (degrees), freq.= ' num2str(fr) ' Hz']);xlabel('r (m)');ylabel('z (m)')
axis equal;axis([0 max(rzb(:,1)) min(rzb(:,2)) max(rzb(:,2))])
colorbar

figure % 3D mesh
surf(RRb,ZZb,abs(p_field))
%shading interp

