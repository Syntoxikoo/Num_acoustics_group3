% Calculates scattering by a half-sphere on a reflecting plane

% Plane wave coming from the +X direction
% The pressure on points over a close arc is calculated

clear

% Geometrical parameters:
R=1;                    % Radius of the sphere
A0=1;                   % Wave amplitude

% AMBIENT CONDITIONS
pa = 101325;          % Atmosferic pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% General parameters:
fr=100;                % Frequency (Hz)
kp=2*pi*fr/c;          % Wavenumber (1/m)
Rfp=R*1.001;           % Radius of the arc of field points (half a circle in theta)
AngleS=90*(pi/180);    % Angle of incidence of a plane wave in the x-z plane with respect to +z, radians (0 to 90, with a reflecting plane)
PlaneON=1;             % If 1, a reflecting plane at z=0 is defined

% GENERATION OF THE DOMAIN GEOMETRY
% Load mesh
fileGEOMsphere='SphereTest2.geo';
[nodesTMP,elementsTMP,elementsQUAD]=readgeomGMSH([fileGEOMsphere '.msh']); % Load mesh from GMSH file
[nodesb,elementsb,segms]=meshcheck(nodesTMP,elementsTMP,0,1); clear nodesTMP elementsTMP;
% Change radius
nodesb(:,1:3)=nodesb(:,1:3)*R;
M=size(nodesb,1); N=size(elementsb,1);


% CHIEF points inside the interior domain:
xyzb_chief=[linspace(-R*0.3,R*0.6,5)' linspace(-R*0.65,R*0.5,5)'  linspace(R*0.11,R*0.55,5)'  ones(5,1)];
hold on; plot3(xyzb_chief(:,1),xyzb_chief(:,2),xyzb_chief(:,3),'md')

% <<<<<<<<<<<<< Run to this point to see the geometry

% Variables used in connection with the reflecting plane:
[nodesbnd,elementsbnd,xyzdum,topodum,xyznodum]=nodummy3D(nodesb,elementsb,PlaneON);
Mnd=size(nodesbnd,1);


% plane wave in negative x-direction
sources=[NaN AngleS 0 A0 0]; 
pI=incoming(sources,[nodesbnd ; xyzb_chief],kp,PlaneON);


% Calculate the BEM matrices and solve the pressures on the surface
[A,B,CConst]=TriQuadEquat(nodesb,elementsb,kp,1,1e-6,PlaneON);  % triangular and/or quadrilateral
% add CHIEF points
[Aex,Bex,Cex]=point(nodesb,elementsb,kp,xyzb_chief,1,1e-6,PlaneON);
A=[A;Aex];B=[B;Bex]; % Extend A and B matrices to include CHIEF points' coefficients
psc=A\(-4*pi*pI);
% clear A B

% Surface pressure, 3D plot:
plotresult(nodesbnd,elementsbnd,psc,1)
axis equal


% Field points
Mfp=25;  % Number of field points
theta=linspace(0,pi,Mfp)';
xyzFP=Rfp*[cos(theta) sin(theta)*sin(0) sin(theta)*cos(0)];
hold on; plot3(xyzFP(:,1), xyzFP(:,2), xyzFP(:,3),'b.')
% Calculate field points
[Afp,Bfp,Cfp]=point(nodesb,elementsb,kp,xyzFP,1,1e-6,PlaneON); % triangular and/or quadrilateral
pIfp=incoming(sources,xyzFP,kp,PlaneON);
pscFP=(Afp*psc)./Cfp + pIfp;


% Analytical solutions:
pscAN = PlaneWaveScatSphere(kp,R,R,acos(-nodesbnd(:,1)/R),1e-6);
pscFPAN = PlaneWaveScatSphere(kp,R,Rfp,theta(end:-1:1),1e-6);  

plotresult(nodesbnd,elementsbnd,psc,1)


% Plot solution (scattering sphere)
figure;
subplot(2,1,1)
plot(theta*180/pi,abs(pscFPAN),'-ko',theta*180/pi,abs(pscFP)/2,'rx'); grid
xlabel('Angle');ylabel('|pressure on field points|');
legend('Analytical','BEM')
title(['Scattering sphere, k=' num2str(kp)]);
subplot(2,1,2)
plot(acos(nodesbnd(:,1)/R)*180/pi,abs(pscAN),'ko',acos(nodesbnd(:,1)/R)*180/pi,abs(psc)/2,'rx'); grid
xlabel('Angle');ylabel('|pressure on the surface|');
legend('Analytical','BEM')

