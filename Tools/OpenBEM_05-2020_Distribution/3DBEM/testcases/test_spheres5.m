% Calculates radiation and scattering by an sphere:
%   -Radiation of a pulsating sphere
%   -Radiation of a first-order oscillating sphere
%   -Scattering of a plane wave by a sphere

% The pressure on points over a a radial direction is calculated
clear

R=1;               % Radius of the sphere
u0=1;              % Maximum velocity amplitude
k=4;               % Wavenumber, m^(-1)
Tole=1e-6;         % Tolerance for the near-singular check. It must be of the order
                   % of the smallest relative distance
phifp=20*pi/180;   % Phi angle (radians) of the line of field points
thetafp=20*pi/180; % Theta angle (radians) of the line of field points
Rfp=[R*(1+1e-6) R*(1+5e-1)];   % Radial range of the field points

nsingON=1;         % Deal with near-singular integrals

% AMBIENT CONDITIONS 
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;              % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% Spheres
% Read nodes and topology. Sphere of radius 1 m.
%nodes=readnodes('sphere.nod');
%elements=readelements('sphere.ele');
%nodes=readnodes('sph4quad.nod');
%elements=readelements('sph4quad.ele');
%nodes=readnodes('sph8quad.nod');
%elements=readelements('sph8quad.ele');
%nodes=readnodes('sphere4tri.nod');
%elements=readelements('sphere4tri.ele');
% for more nodes use sphere8 or sphere16
% for quadratic elements use sphqua4
[nodes,elements,elementsQUAD]=readgeomGMSH(['sphere2.geo.msh']); % Load mesh from GMSH file


% Change radius
nodes(:,1:3)=nodes(:,1:3)*R;

% check geometry and add body numbers
[nodesb,topologyb,segments]=meshcheck(nodes,elements);

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
Mfp=50;  % Number of field points
RR=linspace(Rfp(1),Rfp(2),Mfp)';
xyzFP=RR*[sin(thetafp)*cos(phifp) sin(thetafp)*sin(phifp) cos(thetafp)];
hold on; plot3(xyzFP(:,1), xyzFP(:,2), xyzFP(:,3),'b.')
% Calculate field points
[Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON,Tole); % triangular and/or quadrilateral
pzFP=(Afp*pz+j*k*rho*c*Bfp*vz)./Cfp;
pfFP=(Afp*pf+j*k*rho*c*Bfp*vf)./Cfp;
pscFP=(Afp*psc)./Cfp + exp(i*k*xyzFP(:,3));


% Analytical solutions:
nn=find(RR>R);
pzAN = SphereZerothOrder(k,c,rho,R,u0,RR(nn)); 
pfAN = SphereFirstOrder(k,c,rho,R,u0,[RR(nn) thetafp*ones(length(nn),1)]);
pscAN = PlaneWaveScatSphere(k,R,RR(nn),thetafp,Tole);


% Plot solution (zeroth-order pulsating sphere)
figure;
plot(RR(nn),abs(pzAN),'-k',RR,abs(pzFP),'rx'); grid
xlabel(['Radial direction, theta= ' num2str(thetafp) ' phi= ' num2str(phifp)]);ylabel('|pressure on field points|');
legend('Analytical','Modified BEM')
title(['Zeroth-order pulsating sphere, k=' num2str(k)]);


% Plot solution (first-order oscillating sphere)
figure;
plot(RR(nn),abs(pfAN),'-k',RR,abs(pfFP),'rx'); grid
xlabel(['Radial direction, theta= ' num2str(thetafp) ' phi= ' num2str(phifp)]);ylabel('|pressure on field points|');
legend('Analytical','Modified BEM')
title(['First-order oscillating sphere, k=' num2str(k)]);


% Plot solution (scattering sphere)
figure;
plot(RR(nn),abs(pscAN),'-k',RR,abs(pscFP),'rx'); grid
%plot(RR(nn),abs(pscAN),'-k',RR,abs(pscFP),'rx',RR,abs(pscFP1),'bo'); grid
%xlabel(['Radial direction']);ylabel('|pressure on field points|');
xlabel(['Radial direction, theta= ' num2str(thetafp) ' phi= ' num2str(phifp)]);ylabel('|pressure on field points|');
%legend('Analytical','Standard BEM','Modified BEM')
legend('Analytical','Modified BEM')
title(['Scattering sphere, k=' num2str(k)]);

