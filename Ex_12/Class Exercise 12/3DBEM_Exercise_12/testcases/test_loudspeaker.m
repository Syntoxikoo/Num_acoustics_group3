% Calculates the sound field around a loudpeaker box

% The pressure on points over an arc is calculated to get the directivity
clear

u0=1;              % Velocity amplitude
f=1000;             % Frequency, Hz
nsingON=1;         % Deal with near-singular integrals
lx = 0.135;ly = 0.2;lz = 0.22;        % Box dimensions, m (rectangular)
Rd = 0.05;         % Radius of the diaphragm, assumed centered on the front face (y-z plane) of te box
Dz = 0.03;         % Depth at the dust cap rim
Rfp=2;%ly*2;          % Radius of the arc of field points (half a circle in theta)

% AMBIENT CONDITIONS
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 
k=2*pi*f/c;        % Wavenumber, m^(-1)

% Read nodes and topology. Sphere of radius 1 m.
[nodes,elements,elementsQUAD]=readgeomGMSH(['LoudspeakerBox.msh']); % Load mesh from GMSH file

% check geometry and add body numbers
[nodesb,topologyb,segments]=meshcheck(nodes,elements,0,1);
hold on
M=size(nodesb,1); N=size(topologyb,1);

% Define CHIEF points (random)
Nch=5; % number of points
CHIEFpoint=[linspace(-Dz,-lx,Nch)' linspace(0,ly,Nch)' linspace(0,lz,Nch)'];

% Normal vectors
[nvect]=normals(nodesb,topologyb,'n');

% Velocity of the diaphragm. Convert velocity in the x direction to the normal component
vz=zeros(M,1);
nn=find(sqrt((nodesb(:,2)-ly/2).^2+(nodesb(:,3)-lz/2).^2)<=Rd+eps & nodesb(:,1)>=-Dz-eps);
vx(nn,1)=u0;
vn=dot([vx zeros(M,1) zeros(M,1)]',nvect')'; 
fc=0.3; quiver3(nodesb(:,1),nodesb(:,2),nodesb(:,3),nvect(:,1).*vn*fc,nvect(:,2).*vn*fc,nvect(:,3).*vn*fc,'-r');


% Calculate the BEM matrices and solve the pressures on the surface
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,nsingON);  % triangular and/or quadrilateral
disp(['Condition numbers, A: ' num2str(cond(A)) '  B: ' num2str(cond(B))])

% Include CHIEF points:
[Aex,Bex,Cex]=point(nodesb,topologyb,k,CHIEFpoint);
A=[A;Aex];B=[B;Bex]; % Extend A and B matrices to include CHIEF points' coefficients

% Solve system
B=i*k*rho*c*B;
ps=A\(-B*vn);


% Field points
Mfp=36*2+1;  % Number of field points
thetafp=pi/2*ones(Mfp,1);
phifp=linspace(0,2*pi,Mfp)';
xyzFP=Rfp*[sin(thetafp).*cos(phifp) sin(thetafp).*sin(phifp) cos(thetafp)] + ones(Mfp,1)*[0 ly/2 lz/2];
%xyzFP=Rfp*[sin(thetafp).*cos(phifp) sin(thetafp).*sin(phifp) cos(thetafp)];
hold on; plot3(xyzFP(:,1), xyzFP(:,2), xyzFP(:,3),'b.');axis equal
% Calculate field points
[Afp,Bfp,Cfp]=point(nodesb,topologyb,k,xyzFP,nsingON); % triangular and/or quadrilateral
pzFP=(Afp*ps+j*k*rho*c*Bfp*vn)./Cfp;


% Plot solution
figure;
plot(phifp*180/pi,db(pzFP,20e-6),'-rx'); grid
xlabel('Angle');ylabel('|pressure on field points|')
title(['Loudspeaker box, k=' num2str(k)]);
