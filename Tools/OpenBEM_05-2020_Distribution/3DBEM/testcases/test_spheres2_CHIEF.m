% sphere
% Gives an example of 3-D BEM calculations

% This example adds CHIEF points to the case in "test_spheres1", in order
% to avoid non-uniqueness instabilities for exterior domains.
clear
% Read nodes and topology
% for quadratic elements use sphqua4
% Read nodes and topology. Sphere of radius 1 m.
nodes=readnodes('sphere.nod');
elements=readelements('sphere.ele');
%nodes=readnodes('sph4quad.nod');
%elements=readelements('sph4quad.ele');
%nodes=readnodes('sph8quad.nod');
%elements=readelements('sph8quad.ele');
%nodes=readnodes('sphere4tri.nod');
%elements=readelements('sphere4tri.ele');
% for more nodes use sphere8 or sphere16
% for quadratic elements use sphqua4


k=1; %wavenumber
R=1; % Radius of the sphere

% Change radius
nodes(:,1:3)=nodes(:,1:3)*R;

% check geometry and add body numbers
%[nodesb,topologyb,toposhrinkb,tim,segmopen]=bodyfind(nodes,elements);
[nodesb,topologyb,segments]=meshcheck(nodes,elements);


% CALCULATION AT WELL-BEHAVED FREQUENCY

% Calculate the BEM matrices
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k); 

% The equation for the velocity potential phi is
% 0 = (A-CConst)phi + B v + 4 pi phi^I
% Time dependance exp(i*omega*t)

% plane wave in negative z-direction
phiI=exp(i*k*nodesb(:,3));
phi_sc4=A\(-4*pi*phiI);
plotresult(nodesb,topologyb,abs(phi_sc4));
ana4 = PlaneWaveScatSphere(k,R,R,acos(nodesb(:,3)/R),1e-12);
res4=ana4'-phi_sc4;
% Normalized calculation error: 
error=norm(res4)/norm(ana4)


% Oscillating sphere - axis of symmetry: z
v=(nodesb(:,3));
phi_osc=A\(-B*v);

plotresult(nodesb,topologyb,abs(phi_osc));
b1=j*k^2*exp(j*k)/(k^2-2*j*k-2);
ana_phi_osc=-b1*v*(1+1/j/k)/k*exp(-j*k);
res_rad4=phi_osc-ana_phi_osc;
% Normalized calculation error: 
error_rad=norm(res_rad4)/norm(ana_phi_osc)


% CALCULTAION WITHOUT/WITH CHIEF POINTS AT THE FREQUENCY OF AN INTERNAL RESONANCE

% Calculation with no CHIEF point:
k=pi; % There is an internal resonance at this k value, therefore a non-uniqueness problem
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k); 
pI=exp(i*k*nodesb(:,3));
phi_sc4=A\(-4*pi*pI);
%plotresult(nodesb,topologyb,abs(phi_sc4));
ana4 = PlaneWaveScatSphere(k,R,R,acos(nodesb(:,3)/R),1e-12);
res4=ana4'-phi_sc4;
% Normalized calculation error: 
err_noCHIEF=norm(res4)/norm(ana4)

% Include CHIEF point:
% calculate CHIEF equation for a point in (0,0,0)
CHIEFpoint=[0 0 0];
[Aex,Bex,Cex]=point(nodesb,topologyb,k,CHIEFpoint);
% add CHIEF point
A=[A;Aex];B=[B;Bex]; % Extend A and B matrices to include CHIEF points' coefficients
pI=[pI; exp(i*k*CHIEFpoint(:,3))];
phi_sc4=A\(-4*pi*pI);
plotresult(nodesb,topologyb,abs(phi_sc4));
ana4 = PlaneWaveScatSphere(k,R,R,acos(nodesb(:,3)/R),1e-12);
res4=ana4'-phi_sc4;
% Normalized calculation error: 
err_CHIEF=norm(res4)/norm(ana4)
