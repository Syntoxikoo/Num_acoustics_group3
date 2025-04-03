clear
close all
clc
addpath(genpath("."))
% ADD PARALLEL COMPUTING FOR FASTER PROCESS
%%

radiate = false; % object have infinite impedance
% PREPROCESSING
rho=1.1992;
c=344;

fr=150; ff=1; % (ff is the frequency count when a loop is used)
k=2*pi*fr/c;
% number of elements per wavelength
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points

% ------------------------------- define geometry of the object -------------
Rc=1; % Radius of the cylinder
segments=[-Rc 0 0 Rc ceil(pi/2*Rc*el_wl) Rc 0;...
           0 Rc Rc 0 ceil(pi/2*Rc*el_wl) Rc 0;...
          Rc 0 0 -Rc ceil(pi/2*Rc*el_wl) Rc 0;...
          0 -Rc -Rc 0 ceil(pi/2*Rc*el_wl) Rc 0]; 
[xyb,topology]=nodegen(segments,'y'); 

% ------------------------- generate mesh of field points -------------------

[XX,YY]=meshgrid(-10:espac:10,-5:espac:5);
xy=[XX(1:end)' YY(1:end)'];

%% CALCULATION


inc_pressure=exp(1j*k*xyb(:,1)); % incident pressure at each point of the x coordinates of the object (considering plane wave with A = 1 of x incidence)

if radiate == false
    % calculate coefficient matrix
    A=bem2d(xyb,topology,k(ff));  
    % solve system
    ps = A\(-2*pi*inc_pressure);
elseif radiate == True
    error("bem2d don't solve for radiating object yet")
    ps = A\(-1j*k*rho*c*B*vp);
end

%% Analytical solution for cylinder
pAna=cylscat(k,Rc,xyb(:,1:2),150);
phi=angle(xyb(:,1)+1j*xyb(:,2))';
pAna2 = compute_cylinder_field(xyb(:,1),xyb(:,2),Rc,0,k,1,1000,'cartesian',0,[0,0],true);


% Analytical solution on the field points
[pAnaFP,pAnaFP_inc,pAnaFP_scat]=cylscat(k,Rc,xy,150);
[pAnaFP2,pAnaFP_inc2,pAnaFP_scat2]=compute_cylinder_field(xy(:,1),xy(:,2),Rc,0,k,1,150,'cartesian',0,[0,0],false);
%% BEM solution for the whole field
% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),xy);  


pIfield=exp(1j*k*xy(:,1)); % obtain incident pressure on the field points

% solve the pressure on the field points
pnF = (Ap)*ps./CConst + pIfield;%eye(size(Ap))*inc_pressure; %+ pIfield;

%% plot the pressure on the surface
figure;
plot(rad2deg(phi),abs(pAna), "LineStyle","-",LineWidth=1); hold on; grid on
plot(rad2deg(phi),abs(pAna2), "LineStyle","--",LineWidth=1); 
plot(rad2deg(phi),abs(ps), "LineStyle","-.",LineWidth=1); hold off
xlim([-180 180]);
%% plot the result on the field points 
figure;
subplot(2,1,1)
contourf(XX,YY,reshape(abs(pAnaFP_scat2),size(XX))); grid on 
zlim([-1 1])
subplot(2,1,2)
r = sqrt(xy(:,1).^2 + xy(:,2).^2);
mask = find(r<=Rc);
pnF(mask) = NaN;
contourf(XX,YY,reshape(abs(pnF),size(XX)));

%%
figure;
subplot(2,1,1)
surf(XX,YY,reshape(abs(pAnaFP2),size(XX))); grid on 
% zlim([-1 1])
subplot(2,1,2)
r = sqrt(xy(:,1).^2 + xy(:,2).^2);
mask = find(r<=Rc);
pnF(mask) = NaN;
surf(XX,YY,reshape(abs(pnF),size(XX)));
%%
% ------------------------------- define geometry of the object -------------
lx = 1;
ly=1; 
nbperside = 1;
segments=[-lx/2, ly/2, lx/2, ly/2, nbperside, 0, 0;...
           lx/2 ly/2 lx/2 -ly/2 nbperside 0 0;...
          lx/2 -ly/2 -lx/2 -ly/2 nbperside 0 0;...
          -lx/2 -ly/2 -lx/2 ly/2 nbperside 0 0]; 
[xyb,topology]=nodegen(segments,'y'); 

%% BEM solution for the whole field

inc_pressure=exp(1j*k*xyb(:,1)); % incident pressure at each point of the x coordinates of the object (considering plane wave with A = 1 of x incidence)


A=bem2d(xyb,topology,k(ff));  
% solve system
ps = A\(-2*pi*inc_pressure);
% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k(ff),xy);  


pIfield=exp(1j*k*xy(:,1)); % obtain incident pressure on the field points

% solve the pressure on the field points
pnF = (Ap)*ps./CConst +eye(size(Ap))*inc_pressure; %+ pIfield;

%%
figure;
mask = find(abs(xy) <=lx);
% pnF(mask) =NaN;
surf(XX,YY, reshape(real(pnF),size(XX)))