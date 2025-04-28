% Example of calculation: Rediation from infinite cylinder

% Contains CHIEF points definition

clear

% PREPROCESSING

% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% frequency and wavenumber
fr=300; ff=1; % (ff is the frequency count when a loop is used)
k=2*pi*fr/c;
u0=1;              % Maximum velocity amplitude
anglerad=20;    % Angle of the moving section, degrees

% number of elements per wavelength
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points

% plane admittance (no plane)
sigmaP=NaN;
betaP=betag(sigmaP,fr);
% object admittance
sigmaS=Inf;
betaS=betag(sigmaS,fr)/(rho*c); % de-normalized admittance, as needed to solve the system and domain points


% define geometry of the object
Rc=1; % Radius of the cylinder
segments=[-Rc 0 0 Rc ceil(pi/2*Rc*el_wl) Rc 0;...
           0 Rc Rc 0 ceil(pi/2*Rc*el_wl) Rc 0;...
          Rc 0 0 -Rc ceil(pi/2*Rc*el_wl) Rc 0;...
          0 -Rc -Rc 0 ceil(pi/2*Rc*el_wl) Rc 0]; 
[xyb,topology]=nodegen(segments,'y'); 
M=size(xyb,1);N=size(topology,1);     % M nodes, N elements

% Velocity (section of the cylinder)
nn=find(acos(xyb(:,2)/Rc)<=anglerad*pi/180);
vp=zeros(M,1); vp(nn)=u0;
hold on; plot(xyb(nn,1),xyb(nn,2),'r*')

% CHIEF points:
xyb_chief=[linspace(0,-0.25,10)' linspace(0.5,-0.25,10)' -ones(10,1)];
plot(xyb_chief(:,1),xyb_chief(:,2),'md')


% CALCULATION

% Calculate coefficient matrices
[A,B]=bem2d(xyb,topology,k(ff),betaP,xyb_chief);
B=i*k*rho*c*B;
% solve system
ps=A\(-B*vp); 

% calculation of the pressure on field points:

% generate mesh of field points
[XX,YY]=meshgrid(-3:espac:3,-3:espac:3);
xy=[XX(1:end)' YY(1:end)'];

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),betaP,xy);
% solve the pressure on the field points
pfield=(Ap*ps+j*k*rho*c*Bp*vp)./CConst;

% % Memory saving field point calculation for cases with many field points (replaces the previous two lines)
% for ii=1:size(xy,1);
%    [Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),betaP,xy(ii,:));
%    pfield(ii)=(Ap*ps+j*k*rho*c*Bp*vp)./CConst;
%    disp(['Field point nr. ' num2str(ii) ' of ' num2str(size(xy,1))])
%    clear Ap Bp CConst
% end



% POSTPROCESSING

% plot the pressure on the surface
figure;
plot(1:length(ps),abs(ps),'ko--');
title(['Scattering by a cylinder - Frequency = ' num2str(fr(ff)) ' Hz']);
xlabel('Nodes on the surface');ylabel('Pressure modulus (Pa)')
grid;


% plot the result on the field points
figure;
plotfpoints(XX,YY,abs(pfield));
xlabel('x (m)');
ylabel('y (m)');
title(['Sound Pressure (Pa) - Frequency = ' num2str(fr(ff)) ' Hz']);
axis equal
