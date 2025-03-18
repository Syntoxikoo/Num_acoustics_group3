% Calculates scattering by a sphere above a rectangular sample on a reflecting plane:

% Plane wave coming from the +X direction

clear

% Geometrical parameters:
R=1;                    % Radius of the sphere
Hgt = 2;                % Distance of the centre of the sphere above the origin
Lx=5;                   % Depth of the sample
Ly=4;                   % Width of the sample
Lz=0.5;                 % Height of the sample
A0=1;                   % Wave amplitude

% AMBIENT CONDITIONS
pa = 101325;          % Atmosferic pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% General parameters:
fr=150;                % Frequency (Hz)
kp=2*pi*fr/c;          % Wavenumber (1/m)
AngleS=90*(pi/180);    % Angle of incidence of a plane wave in the x-z plane with respect to +z, radians (0 to 90, with a reflecting plane)
PlaneON=1;             % If 1, a reflecting plane at z=0 is defined

% GENERATION OF THE DOMAIN GEOMETRY
% Load mesh
% fileGEOMsphere='SphereTest4';
fileGEOMsphere='SphereTest3.geo';
[nodesTMP,elementsTMP,elementsQUAD]=readgeomGMSH([fileGEOMsphere '.msh']); % Load mesh from GMSH file
[nodesb,elementsb,segms]=meshcheck(nodesTMP,elementsTMP,0,1); clear nodesTMP elementsTMP;
% Change radius
nodesb(:,1:3)=nodesb(:,1:3)*R;
M=size(nodesb,1); N=size(elementsb,1);
Lz=max(nodesb(nodesb(:,3)<(Hgt-R-eps),3)); % Adjust Lz to different cases

% CHIEF points inside the interior domain:
xyzb_chief=[linspace(-R*0.3,R*0.6,5)' linspace(-R*0.65,R*0.5,5)'  linspace(R*0.11,R*0.55,5)'+Hgt  ones(5,1) ;...
            linspace(-Lx*0.3,Lx*0.4,5)' linspace(-Ly*0.2,Ly*0.35,5)'  linspace(Lz*0.11,Lz*0.9,5)'  ones(5,1)];
hold on; plot3(xyzb_chief(:,1),xyzb_chief(:,2),xyzb_chief(:,3),'md')

% <<<<<<<<<<<<< Run to this point to see the geometry

% Variables used in connection with the reflecting plane:
[nodesbnd,elementsbnd,xyzdum,topodum,xyznodum]=nodummy3D(nodesb,elementsb,PlaneON);
Mnd=size(nodesbnd,1);


% plane wave in negative x-direction
sources=[NaN AngleS 0 A0 0]; 
pI=incoming(sources,[nodesbnd ; xyzb_chief],kp,PlaneON);
%pI=exp(i*kp*nodesbnd(:,1)); % Amplitude = 1


% Admittance of the sample:
Ismp=find(nodesbnd(:,3)<=Lz+eps);
plotsurf(nodesbnd,elementsbnd,[],Ismp)
Y0=1/(rho*c);
Y0=-Y0;  % The sign of the normal vector and the incident particle velocity is opposite
Ys=zeros(Mnd,1); Ys(Ismp)= Y0;

% Calculate the BEM matrices and solve the pressures on the surface
% load TMP A B;
[A,B,CConst]=TriQuadEquat(nodesb,elementsb,kp,1,1e-6,PlaneON);  % triangular and/or quadrilateral
% add CHIEF points
[Aex,Bex,Cex]=point(nodesb,elementsb,kp,xyzb_chief,1,1e-6,PlaneON);
A=[A;Aex];B=[B;Bex]; % Extend A and B matrices to include CHIEF points' coefficients
B=i*kp*rho*c*B;
psc=(A+B*diag(Ys))\(-4*pi*pI);
% clear A B


% Surface pressure, 3D plot:
figure
plotresult(nodesbnd,elementsbnd,psc,1)
title(['Frequency: ' num2str(fr) ' Hz, Impedance ' num2str(1/Y0)])
axis equal

% Field points on a plane in front of the sphere
xi=15;xp=linspace(R,Lx*0.7,xi);
zi=10;zp=linspace(Lz*1.01,Hgt+R,zi);
[XX,ZZ]=meshgrid(xp,zp);
xyzFP=[XX(1:end)' zeros(xi*zi,1) ZZ(1:end)'];
hold on; plot3(xyzFP(:,1), xyzFP(:,2), xyzFP(:,3),'b.')

% Calculate field points
% load TMP Afp Bfp Cfp
[Afp,Bfp,Cfp]=point(nodesb,elementsb,kp,xyzFP,1,1e-6,PlaneON); % triangular and/or quadrilateral
pIfp=incoming(sources,xyzFP,kp,PlaneON);
pscFP=((Afp + i*kp*rho*c*Bfp*diag(Ys))*psc)./Cfp + pIfp;

% Pressure on grid of field points
figure
surf(XX,ZZ,abs(reshape(pscFP,zi,xi)));
axis equal
title(['Frequency: ' num2str(fr) ' Hz, Impedance ' num2str(1/Y0)])
rotate3d on
colorbar
