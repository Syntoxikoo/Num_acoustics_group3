% Example of calculation: Rectangular cavity

clear

% PREPROCESSING

% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% define geometry of the rectangular cavity
lx=1; % Heigth
ly=1.2; % Width

% frequency and wavenumber
m=2;n=3;
fr=c/2*sqrt((m/lx)^2+(n/ly)^2); ff=1; % (ff is the frequency count when a loop is used)
k=2*pi*fr/c;

% number of elements per wavelength
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl/4;       % Spacing of field points

% plane admittance (no plane)
sigmaP=NaN;
betaP=betag(sigmaP,fr);

% object admittance
sigmaS=Inf;
betaS=betag(sigmaS,fr)/(rho*c); % de-normalized admittance, as needed to solve the system and domain points

% define geometry of the object
segments=[0 0 0 ly ceil(ly*el_wl) 0 0;...
          0 ly lx ly ceil(lx*el_wl) 0 0;...
          lx ly lx 0 ceil(ly*el_wl) 0 0;...
          lx 0 0 0 ceil(lx*el_wl) 0 0]; 
[xyb,topology]=nodegen(segments,'y'); 
Ynodes=betaS*ones(size(xyb,1),1);

% Interior problem
xyb(:,end)=-xyb(:,end);
topology(:,end)=-topology(:,end);


% CALCULATION


% line source position
position=[0.08 0.08];

% obtain incident pressure
[G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k(ff),position,betaP(ff),xyb(:,1),xyb(:,2));
inc_pressure=i/4*(G0dir+G0ref)+Pbeta;

% calculate coefficient matrix
[A,B]=bem2d(xyb,topology,k(ff),betaP);
% cond(A)

% solve system
ps=(A+j*k*rho*c*B*diag(Ynodes))\(-2*pi*inc_pressure);

% calculation of the pressure on field points:

% generate mesh of field points
[XX,YY]=meshgrid(linspace(lx/1000,lx-lx/1000,round(lx/espac)),linspace(ly/1000,ly-ly/1000,round(ly/espac)));
xy=[XX(1:end)' YY(1:end)'];

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),betaP,xy);

% obtain incident pressure on the field points
[G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k(ff),position,betaP(ff),xy(:,1),xy(:,2));
pIfield=i/4*(G0dir+G0ref)+Pbeta;
% solve the pressure on the field points
pfield=((Ap+j*k*rho*c*Bp*diag(Ynodes))*ps./CConst+pIfield).';


% POSTPROCESSING

% plot the pressure on the surface
figure;
plot(1:length(ps),abs(ps),'ko--')
title(['Sound pressure on the boundary - Frequency = ' num2str(fr(ff)) ' Hz']);
xlabel('Surface nodes');ylabel('Pressure modulus (Pa)')
grid;


% plot the result on the field points
figure;
plotfpoints(XX,YY,abs(pfield));
xlabel('x (m)');
ylabel('y (m)');
title(['Sund pressure (Pa) - Frequency = ' num2str(fr(ff)) ' Hz']);
axis equal
