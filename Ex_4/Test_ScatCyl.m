% Example of calculation: Infinite cylinder in free field

clear

% PREPROCESSING

% Constants
rho=1.1992;
c=344;

% frequency and wavenumber
fr=150; ff=1; % (ff is the frequency count when a loop is used)
k=2*pi*fr/c;

% number of elements per wavelength
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points


% define geometry of the object
Rc=1; % Radius of the cylinder
segments=[-Rc 0 0 Rc ceil(pi/2*Rc*el_wl) Rc 0;...
           0 Rc Rc 0 ceil(pi/2*Rc*el_wl) Rc 0;...
          Rc 0 0 -Rc ceil(pi/2*Rc*el_wl) Rc 0;...
          0 -Rc -Rc 0 ceil(pi/2*Rc*el_wl) Rc 0]; 
[xyb,topology]=nodegen(segments,'y'); 


%% CALCULATION

% obtain incident pressure
inc_pressure=exp(1j*k*xyb(:,1));

% calculate coefficient matrix
A=bem2d(xyb,topology,k(ff));  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (1) How to calculate the pressure on the surface?

% Analytical solution
pAna=cylscat(k,Rc,xyb(:,1:2),150);
phi=angle(xyb(:,1)+1j*xyb(:,2))';


% Calculation of the pressure on field points:

% generate mesh of field points
[XX,YY]=meshgrid(-10:espac:10,-5:espac:5);
xy=[XX(1:end)' YY(1:end)'];

% obtain incident pressure on the field points
pIfield=exp(1j*k*xy(:,1));

% Analytical solution on the field points
[pAnaFP,pAnaFP_inc,pAnaFP_scat]=cylscat(k,Rc,xy,150);

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),xy);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the pressure on the field points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (2) Use the matrices and boundary results to get the pressure on field points

ps = A\(-2*pi*inc_pressure);
% POSTPROCESSING
pnF = (Ap)*ps./CConst + pIfield;%eye(size(Ap))*inc_pressure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (3) Show the result and compare with analytical data
%%
% plot the pressure on the surface
figure;
plot(rad2deg(phi),abs(pAna), "LineStyle","-",LineWidth=1); hold on; grid on
plot(rad2deg(phi),abs(ps), "LineStyle","--",LineWidth=1); hold off
xlim([-180 180]);
% plot the result on the field points 
figure;
subplot(2,1,1)
surf(XX,YY,reshape(real(pAnaFP_scat),size(XX))); grid on 
subplot(2,1,2)
surf(XX,YY,reshape(real(pnF),size(XX)));