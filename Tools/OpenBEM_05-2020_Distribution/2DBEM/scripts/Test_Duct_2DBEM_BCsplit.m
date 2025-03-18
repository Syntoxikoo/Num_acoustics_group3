% Example of calculation: Rectangular box of plane wave with mixed B.C 

% BC splitting: corner nodes have splitted boudary conditions. The B matrix
% has extra columns to allow it. See PMJ thesis p. 65, sec. 4.4.5.

clear

% PREPROCESSING

% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature 
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% define geometry of the rectangular cavity
ly=0.05; % Heigth
lx=0.3; % Width

% frequency and wavenumber
fr=200;
k=2*pi*fr/c;

% number of elements per wavelength
el_wl=10*max(fr)/c;   % Minimum mesh density as a function of the frequency
espac=1/el_wl/20;       % Spacing of field points
 
% plane admittance (no plane)
sigmaP=NaN;
betaP=betag(sigmaP,fr);

% define geometry of the object
segments=[0 0 0 ly 5 0 0;...
          0 ly lx ly 30 0 0;...
          lx ly lx 0 5 0 0;...
          lx 0 0 0 30 0 0]; 
[xyb,topology]=nodegen(segments,'y'); 

% Interior problem
xyb(:,end)=-xyb(:,end);
topology(:,end)=-topology(:,end);
M=size(xyb,1);

% BC splitting variables, see help in bound2D
BCelem=[ones(5,1) ; ones(30,1)*2 ; ones(5,1)*3 ; ones(30,1)*2]; % assign the elements to one of the three BC regions: lids and wall.
[BCtopo,BCnodeA,BCnodeB]=bound2D(xyb,topology,BCelem,'y');
Ms=size(BCnodeB,1); % number of expanded columns in the B matrix


% Impedance at the lid(s);
IItmp=zeros(Ms,1);
IItmp(BCnodeB(:,2)==1,1)=-1/(rho*c); % rho*c impedance at the emitting end (x=0)
IItmp(BCnodeB(:,2)==3,1)=-1/(rho*c); % rho*c impedance at the receiving end (x=lx)
II=zeros(Ms,M); % Arrange admittances in a non square matrix
for yy=1:Ms % Create row-expanded admittance matrix
    II(yy,BCnodeB(yy,1))=IItmp(yy);
end
clear IItmp


% Piston velocity. Used as excitation.
Vn=zeros(Ms,1); % Combines with the expanded B matrix 
Vn(BCnodeB(:,2)==1,1)=2/(rho*c); % Velocity amplitude corresponding to an incident plane wave of amplitude 1

% CALCULATION


% obtain incident pressure
% inc_pressure=exp(1i*k.*xyb(:,1));

% calculate coefficient matrix (topology matrix is modified to indicate BC splitting, see bem2d help)
[A,B]=bem2d(xyb,[topology(:,1:end-1) topology(:,end)+j*BCelem],k,betaP);
cond(A)

% solve system (rewritten)
ps=(A+1j*k*rho*c*B*II)\(-1j*k*rho*c*B*Vn); % Note that B and II are not square, and Vn matches th expanded coluns of B
% clear Vn
% Vn=II*ps./rho/c;



% calculation of the pressure on field points:

% generate mesh of field points
[XX,YY]=meshgrid(espac:espac:lx,espac:espac:ly);
xy=[XX(1:end)' YY(1:end)'];

% calculate corresponding rows of coefficients (topology matrix is modified to indicate BC splitting, see fieldpoints help)
[Ap,Bp,CConst]=fieldpoints(xyb,[topology(:,1:end-1) topology(:,end)+j*BCelem],k,betaP,xy);

% solve the pressure on the field points
pfield=(Ap*ps+1i*k*rho*c*Bp*(II*ps+Vn))./CConst;


% POSTPROCESSING



% plot the pressure on the surface
figure;
subplot(3,1,1)
plot(1:length(ps),abs(ps),'ko--')
title(['Sound pressure on the boundary - Frequency = ' num2str(fr) ' Hz']);
xlabel('Surface nodes');ylabel('Pressure modulus (Pa)')
grid;
subplot(3,1,2)
plot(1:length(Vn),abs(Vn*rho*c),'ko--')
title(['Normal velocity on the boundary - Frequency = ' num2str(fr) ' Hz']);
xlabel('Surface nodes');ylabel('Velocity (m/s)')
grid;
subplot(3,1,3)
Vn2=zeros(M,1);
Vn2(BCnodeB(:,1),1)=Vn+II*ps;
plot(1:length(ps),abs(Vn2./ps*rho*c),'ko--')
title(['Boundary admittance - Frequency = ' num2str(fr) ' Hz']);
xlabel('Surface nodes');ylabel('Admitttance (V/p)')
grid;


% plot the result on the field points
figure;
plotfpoints(XX,YY,abs(pfield));
xlabel('x (m)');
ylabel('y (m)');
title(['Sound pressure (Pa) - Frequency = ' num2str(fr) ' Hz']);
axis equal



% Animated pressure over length
figure;  % screen animation
phase=0;
while 1 % Execution may be stopped using Ctr-C
    clf;%figure(hh)
    phase=phase+pi/180*10;    if phase>=2*pi-eps, phase=0;end;
    plotfpoints(XX,YY,real(pfield*exp(j*phase)));
    caxis([-max(max(abs(pfield))) max(max(abs(pfield)))])
    xlabel('x (m)');
    ylabel('y (m)');
    title(['Sound pressure (Pa) - Frequency = ' num2str(fr) ' Hz']);
    axis equal
    drawnow;  %   pause(1)
end
close(hh)




