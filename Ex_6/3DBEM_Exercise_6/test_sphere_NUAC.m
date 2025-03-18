
% Calculates scattering by an sphere of a plane wave

% The pressure on points over a close arc is calculated
clear

R=1;               % Radius of the sphere
k=4;               % Wavenumber, m^(-1)
Rfp=R*1.0001;          % Radius of the arc of field points (half a circle in theta)
Tole=1e-6;         % Tolerance for the near-singular check. It must be of the order
                   % of the smallest relative distance
ToleAna=1e-6;      % Tolerance for the analytical solutions
phifp=20*pi/180;   % Phi angle (radians) of the arc of field points (half a circle in theta)
nsingON=1;         % Deal with near-singular integrals

% AMBIENT CONDITIONS
pa = 101325;         % Atmosferic pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% Sphere
[nodes,elements,elementsQUAD]=readgeomGMSH(['Meta_Duct_2.msh']); % Load mesh from GMSH file


% Change radius
nodes(:,1:3)=nodes(:,1:3)*R;

% check geometry and add body numbers
[nodesb,topologyb,segments]=meshcheck(nodes,elements,0,1);

M=size(nodesb,1); N=size(topologyb,1);

% plane wave in negative z-direction
pI=exp(i*k*nodesb(:,3)); % Amplitude = 1

% Calculate the BEM matrices and solve the pressures on the surface
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON,Tole);  % triangular and/or quadrilateral
disp(['Condition numbers, A: ' num2str(cond(A)) '  B: ' num2str(cond(B))])
psc=A\(-4*pi*pI);
clear A B


% Field points
Mfp=25;  % Number of field points
theta=linspace(0,pi,Mfp)';
xyzFP=Rfp*[sin(theta)*cos(phifp) sin(theta)*sin(phifp) cos(theta)];
hold on; plot3(xyzFP(:,1), xyzFP(:,2), xyzFP(:,3),'b.')
% Calculate field points
[Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON,Tole); % triangular and/or quadrilateral
pscFP=(Afp*psc)./Cfp + exp(i*k*xyzFP(:,3));


% Analytical solutions:
pscAN = PlaneWaveScatSphere(k,R,R,acos(nodesb(:,3)/R),ToleAna);
pscFPAN = PlaneWaveScatSphere(k,R,Rfp,theta,ToleAna);  


% Plot solution (scattering sphere)
figure;
subplot(2,1,1)
plot(theta*180/pi,abs(pscFPAN),'-ko',theta*180/pi,abs(pscFP),'rx'); grid
%plot(theta*180/pi,abs(pscFPAN),'-ko',theta*180/pi,abs(pscFP),'rx',theta*180/pi,abs(pscFP1),'bo'); grid
xlabel('Angle');ylabel('|pressure on field points|');
%legend('Analytical','Standard BEM','Modified BEM')
legend('Analytical','BEM')
title(['Scattering sphere, k=' num2str(k)]);
subplot(2,1,2)
plot(acos(nodesb(:,3)/R)*180/pi,abs(pscAN),'ko',acos(nodesb(:,3)/R)*180/pi,abs(psc),'rx'); grid
xlabel('Angle');ylabel('|pressure on the surface|');
legend('Analytical','BEM')



%%
nsingIP = nsing2dQUAD(4,[0.1 0.1 0.0001], [1 0 0 1;1 1 0 1; 0 1 0 1; 0 0 0 1],1e-6,1);
nsingIP2= nsing2dTRIA(4,[0.1 0.1 0.0001], [1 0 0 1; 0 1 0 1; 0 0 0 1],1e-6,1);

%%
