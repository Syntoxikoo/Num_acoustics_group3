% Example of calculation: Infinite cylinder in free field

clear

% PREPROCESSING

% Constants
rho=1.1992;
c=344;

% frequency and wavenumber

k = linspace(5.52,5.5202,100);
C = zeros(length(k),1);
for ii= 1: length(k)
    fr = k(ii)*c/(2*pi);
    % fr=150; ff=1; % (ff is the frequency count when a loop is used)
    
    
    % number of elements per wavelength
    el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
    espac=1/el_wl;       % Spacing of field points
    
    
    % define geometry of the object
    Rc=1; % Radius of the cylinder
    segments=[-Rc 0 0 Rc ceil(pi/2*Rc*el_wl) Rc 0;...
               0 Rc Rc 0 ceil(pi/2*Rc*el_wl) Rc 0;...
              Rc 0 0 -Rc ceil(pi/2*Rc*el_wl) Rc 0;...
              0 -Rc -Rc 0 ceil(pi/2*Rc*el_wl) Rc 0]; 
    [xyb,topology]=nodegen(segments,'e'); 


    % CALCULATION

    
    % obtain incident pressure
    inc_pressure=exp(1j*k(ii)*xyb(:,1));
    
    % calculate coefficient matrix
    A=bem2d(xyb,topology,k(ii)); 
    C(ii) = cond(A);
    
end
figure;
plot(k,C)
[~,idxk]= max(C);

k = k(idxk);
save('Exercice_5/2DBEM_Exercise_5/k_fail.mat', 'k')
%%
% CALCULATION



fr = k*c/(2*pi);



% number of elements per wavelength
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points


% define geometry of the object
Rc=1; % Radius of the cylinder
segments=[-Rc 0 0 Rc ceil(pi/2*Rc*el_wl) Rc 0;...
           0 Rc Rc 0 ceil(pi/2*Rc*el_wl) Rc 0;...
          Rc 0 0 -Rc ceil(pi/2*Rc*el_wl) Rc 0;...
          0 -Rc -Rc 0 ceil(pi/2*Rc*el_wl) Rc 0]; 
[xyb,topology]=nodegen(segments,'e'); 
% obtain incident pressure
inc_pressure=exp(1j*k*xyb(:,1));

% calculate coefficient matrix
A=bem2d(xyb,topology,k); 

% solve system
ps=A\(-2*pi*inc_pressure);

% Analytical solution
pAna=cylscat(k,Rc,xyb(:,1:2),150);


% calculation of the pressure on field points:

% generate mesh of field points
[XX,YY]=meshgrid(-10:espac:10,-5:espac:5);
xy=[XX(1:end)' YY(1:end)'];

% obtain incident pressure on the field points
pIfield=exp(1j*k*xy(:,1));

% Analytical solution
[pAnaFP,pAnaFP_inc,pAnaFP_scat]=cylscat(k,Rc,xy,150);

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k,xy);
% solve the pressure on the field points
pfield=(Ap*ps./CConst+pIfield).';

% POSTPROCESSING

% plot the pressure on the surface
figure;
plot(1:length(ps),abs(ps),'ko--',1:length(pAna),abs(pAna).','kx-');
title(['Scattering by a cylinder - Frequency = ' num2str(fr(ff)) ' Hz']);
xlabel('Nodes on the surface');ylabel('Pressure modulus (Pa)')
legend('BEM','Analytical')
grid;
% 
% 
% % plot the result on the field points, BEM and analytical
% i_int=find(sqrt(xy(:,1).^2+(xy(:,2).^2))<=Rc); 
% pfield(i_int)=NaN; % remove results inside the cylinder
% pAnaFP(i_int)=NaN; % remove results inside the cylinder
% figure;
% subplot(2,1,1)
% surf(XX,YY,reshape(abs(pfield),size(XX)));
% xlabel('x, m'); ylabel('y, m');zlabel('|p|, Pa');
% title(['Scattering by a cylinder (BEM), ka = ' num2str(k*Rc)])
% axis equal
% subplot(2,1,2)
% surf(XX,YY,reshape(abs(pAnaFP),size(XX)));
% xlabel('x, m'); ylabel('y, m');zlabel('|p|, Pa');
% title(['Scattering by a cylinder (Analytical), ka = ' num2str(k*Rc)])
% axis equal
% rotate3d on;
% 
% % 
% % % ANIMATED results, pressure as z value
% % figure
% % hh=figure;  % screen animation
% % ph=0; AziEle=[-40,30]; MaxMin=[-2 +2];
% % while 1 % Execution may be stopped by closing the figure or using Ctr-C
% %     figure(hh)
% %     ph=ph+pi/180*10;    if ph>=2*pi-eps, ph=0;end; % sweep over angles
% %     subplot(2,2,1)
% %     surf(XX,YY,reshape(real(pfield*exp(-j*ph)),size(XX)));
% %     xlabel('x, m'); ylabel('y, m');zlabel('Re[p_{total}], Pa');
% %     title(['Cylinder (BEM), ka = ' num2str(k*Rc)])
% %     axis equal
% %     axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2)) MaxMin])
% %     rotate3d on; view(AziEle)
% %     subplot(2,2,2)
% %     surf(XX,YY,reshape(real((pfield-pIfield.')*exp(-j*ph)),size(XX)));
% %     xlabel('x, m'); ylabel('y, m');zlabel('Re[p_{scattered}], Pa');
% %     title(['Cylinder (BEM), ka = ' num2str(k*Rc)])
% %     axis equal
% %     axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2)) MaxMin])
% %     rotate3d on; view(AziEle)
% %     subplot(2,2,3)
% %     surf(XX,YY,reshape(real(pAnaFP*exp(-j*ph)),size(XX)));
% %     xlabel('x, m'); ylabel('y, m');zlabel('Re[p_{total}], Pa');
% %     title(['Cylinder (Analytical), ka = ' num2str(k*Rc)])
% %     axis equal
% %     axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2)) MaxMin])
% %     rotate3d on; view(AziEle)
% %     subplot(2,2,4)
% %     surf(XX,YY,reshape(real(pAnaFP_scat*exp(-j*ph)),size(XX)));
% %     xlabel('x, m'); ylabel('y, m');zlabel('Re[p_{scattered}], Pa');
% %     title(['Cylinder (Analytical), ka = ' num2str(k*Rc)])
% %     axis equal
% %     axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2)) MaxMin])
% %     rotate3d on; view(AziEle)
% %     drawnow;  %   pause(1)
% % end
% % 
% 
%% With chief point

load("Exercice_5/2DBEM_Exercise_5/k_fail.mat")
c = 343;
fr = k*c/(2*pi);


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
% obtain incident pressure
ChiefPoint = [0.5,-0.2];
inc_pressure=exp(1j*k*xyb(:,1));
pIch=exp(1j*k*ChiefPoint(1));
% calculate coefficient matrix
A=bem2d(xyb,topology,k); 

[Apchief,~,~]=fieldpoints(xyb,topology,k,ChiefPoint);
newA = [A;Apchief];
newP = [inc_pressure;pIch];

%%
% solve system
ps=[A;Apchief]\(-2*pi*[inc_pressure;pIch]);

% Analytical solution
pAna=cylscat(k,Rc,xyb(:,1:2),150);


% calculation of the pressure on field points:

% generate mesh of field points
[XX,YY]=meshgrid(-10:espac:10,-5:espac:5);
xy=[XX(1:end)' YY(1:end)'];

% obtain incident pressure on the field points
pIfield=exp(1j*k*xy(:,1));

% Analytical solution
[pAnaFP,pAnaFP_inc,pAnaFP_scat]=cylscat(k,Rc,xy,150);

% calculate corresponding rows of coefficients
[Ap,Bp,CConst]=fieldpoints(xyb,topology,k,xy);

% solve the pressure on the field points
pfield=(Ap*ps./CConst+pIfield).';

% POSTPROCESSING

% plot the pressure on the surface
figure;
plot(1:length(ps),abs(ps),'ko--',1:length(pAna),abs(pAna).','kx-');
title(['Scattering by a cylinder - Frequency = ' num2str(fr) ' Hz']);
xlabel('Nodes on the surface');ylabel('Pressure modulus (Pa)')
legend('BEM','Analytical')
grid;


%%
save('Exercice_5/2DBEM_Exercise_5/ps_Chief.m', 'ps')