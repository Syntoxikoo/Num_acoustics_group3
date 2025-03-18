
% AxiBEM Example: Sound field in a cylindrical cavity with calculation of a mesh of field points
% Compare with "FEM_CylCavity_Axi_NuAc.m" in axisymmetrical FEM, Matlab PDE toolbox, (Numerical Acoustics) 
% 
% Vicente Cutanda Henriquez 10-2018

clear
% close all

% INPUT DATA

% Cavity dimensions (m):
Rc=0.4; % Cavity Radius
Hc=1; % Cavity Heigth
RTransd=Rc;%0.0175;  % Radius of the vibrating piston (excitation)

% physical constants (MKS units) of air in normal conditions:
% pa = 101325;          % Static pressure (Pa)
% t = 20;               % Temperature (ºC)
% Hr = 50;              % Relative humidity (%)
% [rho,c]=amb2prop(pa,t,Hr,1000); % c: speed of sound (m/s); rho: air density (Kg/m3)
rho=1.1992;
c=344;

% Other parameters:
espac=Rc/10+eps;  % Field point mesh spacing
vampl=1/(rho*c);          % Amplitude of the piston velocity (m/s) 
m=0;              % Axisymmetrical excitation, circumferential mode m=0

% Frequency and wavenumber 
% Reference: "Vibration and Sound", Morse, p.398
Alpha0n = [0 1.2197 2.2331 3.2383 4.2411]; % Only radial and axial waves, all axisymmetrical
nz=2; Alpha=Alpha0n(4);            % Select eigenmode
fr=c/2*sqrt((nz/Hc)^2+(Alpha/Rc)^2); % Frequency (Hz) - Use an analytical resonance
kp=2*pi*fr/c;     % Wavelength (1/m)


% DOMAIN GEOMETRY
% Minimum mesh density as a function of the highest frequency (see nodegen):
el_wl=6*max(fr)/c;   % Elements per meter, six elements per wavelength

% cavity generator
% The objects need to be defined in the clockwise direction by convention
segments=[0 Hc Rc Hc 10 0 el_wl; Rc Hc Rc 0 10 0 el_wl; Rc 0 0 0 10 0 el_wl];


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
disp(['Calculating f= ' num2str(fr) ' Hz'])
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
rr=0:espac/2:max(rzb(:,1)); Nrr=length(rr);
zz=min(rzb(:,2)):espac/2:max(rzb(:,2)); Nzz=length(zz); 
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
title(['Cylinder (BEM), ka = ' num2str(kp*Rc)])
xlabel('r (m)');ylabel('z (m)')
%shading interp

% ANIMATED results, displacement as z value
figure
hh=figure;  % screen animation
ph=0; AziEle=[10,80]; MaxMin=[-max(abs(p_field(1:end))) max(abs(p_field(1:end)))];
while 1 % Execution may be stopped by closing the figure or using Ctr-C
    figure(hh)
    ph=ph+pi/180*10;    if ph>=2*pi-eps, ph=0;end; % sweep over angles
    surf(RR,ZZ,real(p_field*exp(j*ph)));
    xlabel('r, m'); ylabel('z, m');zlabel('Re[p_{total}], Pa');
    title(['Cylinder (BEM), ka = ' num2str(kp*Rc)])
%     axis equal
    axx=axis; axis([axx(1:4) MaxMin])
    rotate3d on; view(AziEle)
    drawnow;  %   pause(1)
end

