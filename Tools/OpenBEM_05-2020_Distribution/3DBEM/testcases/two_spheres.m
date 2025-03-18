% Example of 3-D BEM calculations wth 2 objects (spheres)
clear
% Read nodes and topology
[nodes,elements,elementsQUAD]=readgeomGMSH(['sphere2.geo.msh']); % Load mesh from GMSH file

elements=[elements ; elements+size(nodes,1)];
nodes=[nodes(:,1:2) nodes(:,3)+1.5 ; nodes(:,1:2) nodes(:,3)-1.5];

% check geometry and add body numbers
%[nodesb,topologyb,toposhrinkb,tim,segmopen]=bodyfind(nodes,elements);
[nodesb,topologyb,segments]=meshcheck(nodes,elements,0,1);

k=1; %wavenumber

% Calculate the BEM matrices
[A,B,CConst]=TriQuadEquat(nodesb,topologyb,k,1,1e-6); 

% The equation for the velocity potential phi is
% 0 = (A-CConst)phi + B v + 4 pi phi^I
% Time dependance exp(i*omega*t)

% plane wave in negative z-direction
pI=exp(i*k*nodesb(:,3));
phi_sc4=A\(-4*pi*pI);
plotresult(nodesb,topologyb,abs(phi_sc4));
