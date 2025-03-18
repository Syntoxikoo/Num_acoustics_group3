% Calculates radiation and scattering by an sphere:
%   -Radiation of a pulsating sphere
%   -Radiation of a first-order oscillating sphere
%   -Scattering of a plane wave by a sphere

% The pressure on points over a close arc is calculated
clear

R=1;               % Radius of the sphere
u0=1;              % Maximum velocity amplitude
k=4;               % Wavenumber, m^(-1)
Rfp=R*1.001;          % Radius of the arc of field points (half a circle in theta)
Tole=1e-6;         % Tolerance for the near-singular check. It must be of the order
                   % of the smallest relative distance
phifp=20*pi/180;   % Phi angle (radians) of the arc of field points (half a circle in theta)
nsingON=1;         % Deal with near-singular integrals
% CONDICIONES AMBIENTALES
pa = 101325;         % Presion Atmosferica (Pa)
t = 20;              % Temperatura (ºC)
Hr = 50;             % Humedad relativa (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% Spheres
% Read nodes and topology. Sphere of radius 1 m.
%nodes=readnodes('sphere.nod');
%elements=readelements('sphere.ele');
% nodes=readnodes('sphere8.nod');
% elements=readelements('sphere8.ele');
% nodes=readnodes('sphere16.nod');
% elements=readelements('sphere16.ele');
%  
% nodes=readnodes('sph4quad.nod');
% elements=readelements('sph4quad.ele');
% nodes=readnodes('sph8quad.nod');
% elements=readelements('sph8quad.ele');
% nodes=readnodes('sphere4tri.nod');
% elements=readelements('sphere4tri.ele');
% for more nodes use sphere8 or sphere16
% for quadratic elements use sphqua4
[nodes,elements,elementsQUAD]=readgeomGMSH(['sphere2.geo.msh']); % Load mesh from GMSH file


% Change radius
nodes(:,1:3)=nodes(:,1:3)*R;

% check geometry and add body numbers
[nodesb,topologyb,segments]=meshcheck(nodes,elements,0,1);

M=size(nodesb,1); N=size(topologyb,1);

% Velocity (zeroth-order pulsating sphere)
vz=ones(M,1)*u0;

% Velocity (first-order oscillating sphere)
vf=u0*nodesb(:,3)/R;

% plane wave in negative z-direction
pI=exp(i*k*nodesb(:,3)); % Amplitude = 1

% Calculate the BEM matrices and solve the pressures on the surface
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON,Tole);  % triangular and/or quadrilateral
disp(['Condition numbers, A: ' num2str(cond(A)) '  B: ' num2str(cond(B))])
B=i*k*rho*c*B;
pz=A\(-B*vz); % the three test cases
pf=A\(-B*vf);
psc=A\(-4*pi*pI);
clear A B


% Field points
Mfp=25;  % Number of field points
theta=linspace(0,pi,Mfp)';
xyzFP=Rfp*[sin(theta)*cos(phifp) sin(theta)*sin(phifp) cos(theta)];
hold on; plot3(xyzFP(:,1), xyzFP(:,2), xyzFP(:,3),'b.')
% Calculate field points
[Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON,Tole); % triangular and/or quadrilateral
pzFP=(Afp*pz+j*k*rho*c*Bfp*vz)./Cfp;
pfFP=(Afp*pf+j*k*rho*c*Bfp*vf)./Cfp;
pscFP=(Afp*psc)./Cfp + exp(i*k*xyzFP(:,3));


% Analytical solutions:
pzAN = SphereZerothOrder(k,c,rho,R,u0,R*ones(M,1)); 
pzFPAN = SphereZerothOrder(k,c,rho,R,u0,Rfp*ones(Mfp,1));
pfAN = SphereFirstOrder(k,c,rho,R,u0,[R*ones(M,1) acos(nodesb(:,3)/R)]);
pfFPAN = SphereFirstOrder(k,c,rho,R,u0,[Rfp*ones(Mfp,1) theta]);
pscAN = PlaneWaveScatSphere(k,R,R,acos(nodesb(:,3)/R),Tole);
pscFPAN = PlaneWaveScatSphere(k,R,Rfp,theta,Tole);  


% Plot solution (zeroth-order pulsating sphere)
figure;
subplot(2,1,1)
plot(theta*180/pi,abs(pzFPAN),'-ko',theta*180/pi,abs(pzFP),'-rx'); grid
xlabel('Angle');ylabel('|pressure on field points|');
legend('Analytical','BEM')
title(['Zeroth-order pulsating sphere, k=' num2str(k)]);
subplot(2,1,2)
plot(acos(nodesb(:,3)/R)*180/pi,abs(pzAN),'ko',acos(nodesb(:,3)/R)*180/pi,abs(pz),'rx'); grid
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
plot(acos(nodesb(:,3)/R)*180/pi,abs(pfAN),'ko',acos(nodesb(:,3)/R)*180/pi,abs(pf),'rx'); grid
xlabel('Angle');ylabel('|pressure on the surface|');
legend('Analytical','BEM')


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

