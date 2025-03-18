% sphere
% Gives an example of 3-D BEM calculations
clear
% Read nodes and topology
nodes=readnodes('sphere.nod');
elements=readelements('sphere.ele');
% for more nodes use sphere8 or sphere16
% for quadratic elements use sphqua4

% check geometry and add body numbers
%[nodesb,topologyb,toposhrinkb,tim,segmopen]=bodyfind(nodes,elements);
[nodesb,topologyb,segments]=meshcheck(nodes,elements);

% interior problem
%nodesb(:,end)=-nodesb(:,end);
%elementsb(:,end)=-elementsb(:,end);

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
% ana4=ana_sphere(k,nodesb,1e-12);
% res4=ana4-phi_sc4;
% error=norm(res4)/norm(ana4);


% Oscillating sphere - axis of symmetry: z
v=(nodesb(:,3));
phi_osc=A\(-B*v);

plotresult(nodesb,topologyb,abs(phi_osc));
