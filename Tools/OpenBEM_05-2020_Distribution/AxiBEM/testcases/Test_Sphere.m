% Calculates two radiation examples:
%   -a piston on a sphere
%   -A first-order oscillating sphere
% The pressure at some distant field points is also calculated
clear

R=1;               % Radius of the sphere
u0=1;              % Maximum velocity amplitude
anglepiston=20;    % Piston angle in degrees
Rfp=R*10;          % Radius of the arc of field points (half a circle in the z-y plane)
k=1;               % Wavenumber, m-1
m=0;               % Axisymmetrical excitation, first mode m=0
quadelem=1;        % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)

% AMBIENT CONDITIONS
pa = 101325;          % Atmosferic pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

 

% GENERATION OF THE DOMAIN GEOMETRY
N_nodes=30;
a_p=(90-anglepiston)*pi/180;   % Piston angle in radians
a_tr=(anglepiston/100)*pi/180; % Transition
% segments=[0 R                             R*cos(a_p) R*sin(a_p)            round(N_nodes*anglepiston/180/2) R 0; ...
%           R*cos(a_p) R*sin(a_p)           0 -R                             round(N_nodes*(180-anglepiston)/180/2) R 0]; % No transition
segments=[0 R                             R*cos(a_p) R*sin(a_p)            round(N_nodes*anglepiston/180/2) R 0;
          R*cos(a_p) R*sin(a_p)           R*cos(a_p-a_tr) R*sin(a_p-a_tr)  1 R 0; ...
          R*cos(a_p-a_tr) R*sin(a_p-a_tr) 0 -R                             round(N_nodes*(180-anglepiston)/180/2) R 0]; % With transition element
[rzb,topology]=nodegen(segments,'y',{},quadelem); % nodes and elements
M=size(rzb,1);N=size(topology,1);     % M nodos, N elementos


% % Velocity (piston on sphere)
nn=find(acos(rzb(:,2)/R)<=anglepiston*pi/180);
vp=zeros(M,1); vp(nn)=u0;

% Velocity (first-order oscillating sphere)
vf=u0*rzb(:,2)/R;

hold on; plot(rzb(nn,1),rzb(nn,2),'r*')


% <<<<<<<<<<<<< Run to this point to see the geometry


% BEM MATRICES CALCULATION

[A,B,CConst]=BEMEquat0(rzb,topology,k,m);
disp(['Condition numbers, A: ' num2str(cond(A)) ' and B: ' num2str(cond(B))])
B=i*k*rho*c*B;
% SURFACE SOLUTION
pp=A\(-B*vp); 
pf=A\(-B*vf); % the two test cases


% Field points
Mfp=25;  % Number of field points
theta=linspace(0,pi,Mfp)';
fprzFP=Rfp*[sin(theta) cos(theta)];
% Calculate field points
ppFP=FieldPnt3(fprzFP,pp,vp,rzb,topology,k,m,rho,c);
pfFP=FieldPnt3(fprzFP,pf,vf,rzb,topology,k,m,rho,c);


% Analytical solutions:
pfAN = SphereFirstOrder(k,c,rho,R,u0,[R*ones(M,1) acos(rzb(:,2)/R)]);
pfFPAN = SphereFirstOrder(k,c,rho,R,u0,[Rfp*ones(Mfp,1) theta]);
for ii=1:Mfp
   disp(['Analytical solution (piston), point ' num2str(ii)]);
   ppFPAN(ii) = PistonOnSphere(k*c/(2*pi),c,rho,R,anglepiston*pi/180,u0,Rfp,theta(ii));
end
for ii=1:M
   disp(['Analytical solution (piston), point ' num2str(ii)]);
   ppAN(ii) = PistonOnSphere(k*c/(2*pi),c,rho,R,anglepiston*pi/180,u0,R,acos(rzb(ii,2)/R));
end

% Plot solution (piston on sphere)
figure;
subplot(2,1,1)
plot(theta*180/pi,abs(ppFPAN),'-ko',theta*180/pi,abs(ppFP),'-rx'); grid
xlabel('Angle');ylabel('|pressure on field points|');
legend('Analytical','BEM')
title(['Piston on sphere, k=' num2str(k) ', Rfp=' num2str(Rfp/R) '*R m']);
subplot(2,1,2)
plot(acos(rzb(:,2)/R)*180/pi,abs(ppAN),'ko',acos(rzb(:,2)/R)*180/pi,abs(pp),'rx'); grid
xlabel('Angle');ylabel('|pressure on the surface|');
legend('Analytical','BEM')


% Plot solution (first-order oscillating sphere)
figure;
subplot(2,1,1)
plot(theta*180/pi,abs(pfFPAN),'-ko',theta*180/pi,abs(pfFP),'-rx'); grid
xlabel('Angle');ylabel('|pressure on field points|');
legend('Analytical','BEM')
title(['First-order oscillating sphere, k=' num2str(k)]);
subplot(2,1,2)
plot(acos(rzb(:,2)/R)*180/pi,abs(pfAN),'ko',acos(rzb(:,2)/R)*180/pi,abs(pf),'rx'); grid
xlabel('Angle');ylabel('|pressure on the surface|');
legend('Analytical','BEM')



