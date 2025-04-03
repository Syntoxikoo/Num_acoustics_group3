function [rzb,topology,rzline,IndMb,IndSide,Mmb,M,N,Tm,mu_m,Rad,Hgh,BPrad]=geom_BK4938_9_VER7(MeshDens,SubVer)

% This function produces the geometry and physical parameters of the
% idealized microphone with no back cavity in VCH thesis and KAmpinga
% thesis.

% The input parameter is the mesh density.
% The output contains, for each microphone: node and element matrices,
% indices of the membrane nodes, number of nodes and elements, membrane
% tension and menbrane surface density, in this order.
%
% Rounded edges are used, with adjustable curvature radius.

% Vicente Cutanda Henríquez 04-2011


% Physical dimensions of the microphone's interior, not counting the cavity around.
% Membrane radius is 2 mm, gap thickness is 18e-6 and backplate radius is 1.75 mm
% for microphones 4938 and 4939.
Hgh=18e-6; % Other values: 11e-6 - 24e-6
Rad=2e-3; 
BPrad=2e-3;
Rc=Hgh/5; %5e-5;1e-5; Radius of rounded corners. Hgh=18e-6 a,b: Hgh/5,
ElMem=40;   % elements at the membrane/backplate
ElSide=10;  % elements at the side
ElRound=5;  % elements at the rounded edges
quadelem=1; % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)

% Generator of microphones 4938 and 4839, from B&K drawings.

if SubVer(1)=='a'
    % Serie VER7_a: Uniform mesh
    segments=[0 Hgh Rad Hgh ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='b'
    %Serie VER7_b: Gradual mesh segment close to the rim
    % The number of integration points in IntF1 is increased to 100.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    segments=[0 Hgh Rad*0.8 Hgh ElMem 0 MeshDens;...
        Rad*0.8 Hgh Rad Hgh j*ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 Rad*0.8 0 j*ElMem 0 MeshDens; ...
        Rad*0.8 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='c'
    % Serie VER7_a: Uniform mesh, 3*ElMem
    segments=[0 Hgh Rad Hgh 3*ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 0 0 3*ElMem 0 MeshDens];
elseif SubVer(1)=='d'
    % Serie VER7_a: Uniform mesh, like 'a'.
    % The number of integration points in IntF1 is increased to 100.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    segments=[0 Hgh Rad Hgh ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='e'
    % Serie VER7_a: Gradual mesh, like 'b'.
    % The number of integration points in IntF1 is increased to 200.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    segments=[0 Hgh Rad*0.8 Hgh ElMem 0 MeshDens;...
        Rad*0.8 Hgh Rad Hgh j*ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 Rad*0.8 0 j*ElMem 0 MeshDens; ...
        Rad*0.8 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='f'
    % Serie VER7_a: Gradual mesh, like 'b'.
    % The number of integration points in IntF1 is 200.
    % The order od the gausian integration in IntF1 changes from 10 to 30
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    segments=[0 Hgh Rad*0.8 Hgh ElMem 0 MeshDens;...
        Rad*0.8 Hgh Rad Hgh j*ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 Rad*0.8 0 j*ElMem 0 MeshDens; ...
        Rad*0.8 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='g'
    % Same as "f", but with standard ambient conditions and constants (VTconst)
    segments=[0 Hgh Rad*0.8 Hgh ElMem 0 MeshDens;...
        Rad*0.8 Hgh Rad Hgh j*ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 Rad*0.8 0 j*ElMem 0 MeshDens; ...
        Rad*0.8 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='h'
    % Serie VER7_a: Uniform mesh, like 'a'.
    % The number of integration points in IntF1 is increased to 100.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    % NEW: IntF2 is an earlier version providing better accuracy
    % NEW: IntF1 (nonsingular) is also refined with near-singular integration
    segments=[0 Hgh Rad Hgh ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='i'
    % Serie VER7_a: Gradual mesh, like 'b'.
    % The number of integration points in IntF1 is 200.
    % The order od the gausian integration in IntF1 changes from 10 to 30
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    % NEW: IntF2 is an earlier version providing better accuracy
    % NEW: IntF1 (nonsingular) is also refined with near-singular integration
    segments=[0 Hgh Rad*0.8 Hgh ElMem 0 MeshDens;...
        Rad*0.8 Hgh Rad Hgh j*ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 Rad*0.8 0 j*ElMem 0 MeshDens; ...
        Rad*0.8 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='j'
    % Serie VER7_a: Uniform mesh, like 'a'.
    % The number of integration points in IntF1 is increased to 100.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    % NEW: IntF2 is an earlier version providing better accuracy
    % NEW: IntF1 (nonsingular) is also refined with near-singular integration
    ElMem=20;   % elements at the membrane/backplate
    ElSide=5;  % elements at the side
    ElRound=3;  % elements at the rounded edges
    segments=[0 Hgh Rad Hgh ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='k'
    % Serie VER7_a: Uniform mesh, like 'a'.
    % The number of integration points in IntF1 is increased to 100.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    % NEW: IntF2 is an earlier version providing better accuracy
    % NEW: IntF1 (nonsingular) is also refined with near-singular integration
    ElMem=20;   % elements at the membrane/backplate
    ElSide=5;  % elements at the side
    segments=[0 Hgh Rad Hgh ElMem 0 MeshDens;...
        Rad Hgh Rad 0 ElSide 0 MeshDens;...
        Rad 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='l'
    % Serie VER7_a: Uniform mesh, like 'a'.
    % The number of integration points in IntF1 is increased to 100.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    % NEW: IntF2 is an earlier version providing better accuracy
    % NEW: IntF1 (nonsingular) is also refined with near-singular integration
    ElMem=10;   % elements at the membrane/backplate
    ElSide=5;  % elements at the side
    segments=[0 Hgh Rad Hgh ElMem 0 MeshDens;...
        Rad Hgh Rad 0 ElSide 0 MeshDens;...
        Rad 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='m'
    % Serie VER7_a: Uniform mesh, like 'a'.
    % The number of integration points in IntF1 is increased to 100.
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    % NEW: IntF2 is an earlier version providing better accuracy
    % NEW: IntF1 (nonsingular) is also refined with near-singular integration
    ElMem=10;   % elements at the membrane/backplate
    ElSide=5;  % elements at the side
    ElRound=3;  % elements at the rounded edges
    segments=[0 Hgh Rad Hgh ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 0 0 ElMem 0 MeshDens];
elseif SubVer(1)=='n'
    % Serie VER7_a: Gradual mesh, like 'b' and 'i'. LINEAR elements are used.
    % The number of integration points in IntF1 is 200.
    % The order od the gausian integration in IntF1 changes from 10 to 30
    % The version in 'Version0' is used with modified elliptic integral
    % The bug of changed viscosity coefficients is corrected
    % NEW: IntF2 is an earlier version providing better accuracy
    % NEW: IntF1 (nonsingular) is also refined with near-singular integration
    quadelem=0; % LINEAR elements are used
    segments=[0 Hgh Rad*0.8 Hgh ElMem 0 MeshDens;...
        Rad*0.8 Hgh Rad Hgh j*ElMem 0 MeshDens;...
        Rad Hgh Rad+Rc Hgh-Rc ElRound Rc MeshDens;...
        Rad+Rc Hgh-Rc Rad+Rc 0+Rc ElSide 0 MeshDens;...
        Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
        Rad 0 Rad*0.8 0 j*ElMem 0 MeshDens; ...
        Rad*0.8 0 0 0 ElMem 0 MeshDens];
end


% Serie VER7_c: ??????????????


% % Example with a Bézier corner
% segments=[0 Hgh Rad*0.8 Hgh ElMem 0 MeshDens;...
%           Rad*0.8 Hgh Rad Hgh j*ElMem 0 MeshDens;...
%           Rad+j*Hgh Rad+Rc+j*Hgh Rad+Rc+j*Hgh Rad+Rc+j*(Hgh-Rc) ElRound*10 Rc MeshDens;...
%           Rad+Rc Hgh-Rc Rad+Rc 0+Rc j*ElSide 0 MeshDens;...
%           Rad+Rc 0+Rc Rad 0 ElRound Rc MeshDens;...
%           Rad 0 Rad*0.8 0 j*ElMem 0 MeshDens; ...
%           Rad*0.8 0 0 0 ElMem 0 MeshDens];

      
[rzb,topology,rzline]=nodegen(segments,'y',{},quadelem); 

M=size(rzb,1);N=size(topology,1);

% indicate interior domain
rzb(:,end)=-rzb(:,end);
topology(:,end)=-topology(:,end);


% find nodes on the membrane:
IndMb=find(rzb(:,2)>Hgh-eps & rzb(:,2)<Hgh+eps & rzb(:,1)<=Rad*(1+eps));
Mmb=length(IndMb);
%IndSide=find(rzb(:,1)>(Rad)+eps);
IndSide=find(rzb(:,2)<Hgh-eps & rzb(:,2)>eps);

% specific parameters for the diaphragms, mics [4938 4939], used in solving process.
Tm=3128; %Tm39=1039; %[5700 3450] [4000 1200] tension, originally [3128 1039]. 
% [5700 3450] give vacuum resonances of 60 and 80 kHz, respectively for 4938 and 4839 mics.
mu_m=8300*6.95e-6; %mu_m39=8300*2.35e-6;% surface density, density*thickness 8300*[6.95e-6 2.35e-6]


