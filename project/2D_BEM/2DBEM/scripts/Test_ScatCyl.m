% Example of calculation: Infinite cylinder in free field

clear

% PREPROCESSING

% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% frequency and wavenumber
fr=150; ff=1; % (ff is the frequency count when a loop is used)
k=2*pi*fr/c;

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

% % Interior problem
% xyb(:,end)=-xyb(:,end);
% topology(:,end)=-topology(:,end);

% CHIEF points:
xyb_chief=[linspace(0,-0.25,10)' linspace(0.5,-0.25,10)' -ones(10,1)];
inc_pressure_chief=2*pi*exp(i*k*xyb_chief(:,1));
hold on;plot(xyb_chief(:,1),xyb_chief(:,2),'md')

% CALCULATION

% obtain incident pressure
inc_pressure=2*pi*exp(i*k*xyb(:,1));

% calculate coefficient matrix
A=bem2d(xyb,topology,k(ff),betaP,xyb_chief);

% solve system
% ps=A\(-inc_pressure); % no CHIEF points
ps=A\(-[inc_pressure ; inc_pressure_chief]);

% Analytical solution
divs=(size(xyb,1));
pAna=cylscat(k(ff)*Rc,divs,150,'n');
pAna=pAna(1:divs+1);


% calculation of the pressure on field points:

% generate mesh of field points
[XX,YY]=meshgrid(-10:espac:10,-5:espac:5);
xy=[XX(1:end)' YY(1:end)'];

% obtain incident pressure on the field points
pIfield=exp(i*k*xy(:,1));

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),betaP,xy);
% solve the pressure on the field points
pfield=(Ap*ps./CConst+pIfield).';

% % Memory saving field point calculation for cases with many field points (replaces the previous two lines)
% for ii=1:size(xy,1);
%    [Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),betaP,xy(ii,:));
%    pfield(ii)=(Ap*ps./CConst+pIfield(ii)).';
%    disp(['Field point nr. ' num2str(ii) ' of ' num2str(size(xy,1))])
%    clear Ap CConst
% end



% POSTPROCESSING

% plot the pressure on the surface
figure;
plot(1:length(ps),abs(ps),'ko--',1:length(pAna),abs(pAna).','kx-');
title(['Scattering by a cylinder - Frequency = ' num2str(fr(ff)) ' Hz']);
xlabel('Nodes on the surface');ylabel('Pressure modulus (Pa)')
legend('BEM','Analytical')
grid;


% plot the result on the field points
figure;
plotfpoints(XX,YY,abs(pfield));
xlabel('x (m)');
ylabel('y (m)');
title(['Sound Pressure (Pa) - Frequency = ' num2str(fr(ff)) ' Hz']);
axis equal
