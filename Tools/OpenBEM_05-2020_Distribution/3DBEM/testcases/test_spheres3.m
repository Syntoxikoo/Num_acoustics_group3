% Calculates three examples:
%   -Radiation of a pulsating sphere
%   -Radiation of a first-order oscillating sphere
%   -Radiation of a a piston on a sphere

% The pressure at an arc of field points is also calculated
clear

R=1;               % Radius of the sphere
u0=1;              % Maximum velocity amplitude
anglepiston=20;    % Angle of the piston, degrees
Rfp=R*10;           % Radius of the arc of field points (half a circle in the z-y plane)
k=1;               % Wavenumber, m-1
nsingON=1;         % Deal with near-singular integrals
% CONDICIONES AMBIENTALES
pa = 101325;         % Presion Atmosferica (Pa)
t = 20;              % Temperatura (ºC)
Hr = 50;              % Humedad relativa (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 


% Read nodes and topology. Sphere of radius 1 m.
%nodes=readnodes('sphere.nod');
%elements=readelements('sphere.ele');
%nodes=readnodes('sphere4tri.nod');
%elements=readelements('sphere4tri.ele');
% for more nodes use sphere8 or sphere16
% for quadratic elements use sphqua4
[nodes,elements,elementsQUAD]=readgeomGMSH(['sphere2.geo.msh']); % Load mesh from GMSH file


% Change radius
nodes(:,1:3)=nodes(:,1:3)*R;

% check geometry and add body numbers
%[nodesb,topologyb,toposhrinkb,tim,segmopen]=bodyfind(nodes,elements);
[nodesb,topologyb,segments]=meshcheck(nodes,elements);

M=size(nodesb,1); N=size(topologyb,1);


% Velocity (zeroth-order pulsating sphere)
vz=ones(M,1)*u0;

% Velocity (first-order oscillating sphere)
vf=u0*nodesb(:,3)/R;

% Velocity (piston on sphere)
nn=find(acos(nodesb(:,3)/R)<=anglepiston*pi/180);
vp=zeros(M,1); vp(nn)=u0;

% Calculate the BEM matrices and solve the pressures on the surface
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON); % quadrilateral and/or triangular
B=i*k*rho*c*B;
pz=A\(-B*vz); % the three test cases
pf=A\(-B*vf);
pp=A\(-B*vp); 


% Field points
Mfp=25;  % Number of field points
theta=linspace(0,pi,Mfp)';
xyzFP=Rfp*[zeros(Mfp,1) sin(theta) cos(theta)];
% Calculate field points
[Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON); % quadrilateral
%[Afp,Bfp,Cfp]=point(nodesb,topologyb,topologyb,ones(N,1),k,xyzFP,nsingON);  % triangular
pzFP=(Afp*pz+j*k*rho*c*Bfp*vz)./Cfp;
pfFP=(Afp*pf+j*k*rho*c*Bfp*vf)./Cfp;
ppFP=(Afp*pp+j*k*rho*c*Bfp*vp)./Cfp;


% Analytical solutions:
pzAN = SphereZerothOrder(k,c,rho,R,u0,R*ones(M,1));
pzFPAN = SphereZerothOrder(k,c,rho,R,u0,Rfp*ones(Mfp,1));
pfAN = SphereFirstOrder(k,c,rho,R,u0,[R*ones(M,1) acos(nodesb(:,3)/R)]);
pfFPAN = SphereFirstOrder(k,c,rho,R,u0,[Rfp*ones(Mfp,1) theta]);
for ii=1:Mfp
   disp(['Analytical solution (piston), point ' num2str(ii)]);
   ppFPAN(ii) = PistonOnSphere(k*c/(2*pi),c,rho,R,anglepiston,u0,Rfp,theta(ii));
end
for ii=1:M
   disp(['Analytical solution (piston), point ' num2str(ii)]);
   ppAN(ii) = PistonOnSphere(k*c/(2*pi),c,rho,R,anglepiston,u0,R,acos(nodesb(ii,3)/R));
end

% Plot solution (piston on sphere)
figure;
subplot(2,1,1)
plot(theta*180/pi,abs(ppFPAN),'-ko',theta*180/pi,abs(ppFP)*mean(abs(ppFPAN')./abs(ppFP)),'-rx'); grid
xlabel('Angle');ylabel('|pressure on field points|');
legend('Analytical','BEM')
title(['Piston on sphere, k=' num2str(k)]);
subplot(2,1,2)
plot(acos(nodesb(:,3)/R)*180/pi,abs(ppAN),'ko',acos(nodesb(:,3)/R)*180/pi,abs(pp),'rx'); grid
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

