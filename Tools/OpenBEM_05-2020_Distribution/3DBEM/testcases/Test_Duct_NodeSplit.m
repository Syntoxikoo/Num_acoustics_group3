
% Plane wave propagating in a rectangular duct. The impedance at the two ends
% is set to rho*c to avoid reflections, but can be changed.
%
% One of the lids moves as a piston. The nodes at the rims of the end lids 
% have a split boundary condition,splitting of admittance and velocity.

% VCH, 11-2016

clear, %close all
tic

% Try the same calculation with and without node splitting to see the difference
SplitNodes=1; % If 1, rim nodes at the lids are split; 0, no split

% Frequency
f=1000;

% Dimensions:
lx=20e-2; ly=2e-2; lz=1e-2;

% Mesh file, from GMSH
fileGEOMduct='Meta_Duct_1'; % Rectangular duct with no metamaterial
% fileGEOMduct='Meta_Duct_2'; % Rectangular duct with no metamaterial. Denser mesh

% physical constants of air in normal conditions:
pa = 101325;          % Static pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Relative humidity (%)
[rho,c]=amb2prop(pa,t,Hr,1000);
kp=2*pi*f/c;

% Amplitude of the normal velocity excitation
va=1/(rho*c); % This value should ideally produce a pressure wave of amplitude 1

% Admittance of the anechoic termination. Minus sign due to mismatch with the BEM normal vector.
EndAdm=-1/(rho*c);
ImpE=1; % If set to 1, "EndAdm" impedance is set on both ends of the tube, otherwise only at the receiving end


% DUCT - Create meshes:
[nodes,elementsTRI,elementsQUAD]=readgeomGMSH([fileGEOMduct '.msh']); % Load mesh from GMSH file
[xyzbDC,topologybDC,segms]=meshcheck(nodes,elementsTRI,0,0); clear nodes elementsTRI;
Mdc=size(xyzbDC,1);Ndc=size(topologybDC,1);
% Define interior domain
xyzbDC=[xyzbDC(:,1:end-1) -xyzbDC(:,end)];
topologybDC=[topologybDC(:,1:end-1) -topologybDC(:,end)];

% Resize to prescribed dimensions and set origin
xmax = max(xyzbDC(:,1)); xmin = min(xyzbDC(:,1));
ymax = max(xyzbDC(:,2)); ymin = min(xyzbDC(:,2));
zmax = max(xyzbDC(:,3)); zmin = min(xyzbDC(:,3));
xyzbDC(:,1) = (xyzbDC(:,1) -xmin)*lx/(xmax-xmin); % set origin at one end and scale to the prescribed x dimension
xyzbDC(:,2) = (xyzbDC(:,2) -ymin)*ly/(ymax-ymin); % set origin at one end and scale to the prescribed y dimension
xyzbDC(:,3) = (xyzbDC(:,3) -zmin)*lz/(zmax-zmin); % set origin at one end and scale to the prescribed z dimension



% Find nodes and elements on the excitation side and receiving side (Duct 2)
NodDCe=find(xyzbDC(:,1)<=eps & xyzbDC(:,2)<ly-eps & xyzbDC(:,2)>eps & xyzbDC(:,3)<lz-eps & xyzbDC(:,3)>eps);  Mdce=length(NodDCe);
NodDCr=find(xyzbDC(:,1)>lx-eps & xyzbDC(:,2)<ly-eps & xyzbDC(:,2)>eps & xyzbDC(:,3)<lz-eps & xyzbDC(:,3)>eps);  Mdcr=length(NodDCr);
% NodDCe=find(xyzbDC(:,1)<=eps);  Mdce=length(NodDCe);
% NodDCr=find(xyzbDC(:,1)>lx-eps);  Mdcr=length(NodDCr);


ElDCe=find( xyzbDC(topologybDC(:,1),1)<=eps & xyzbDC(topologybDC(:,2),1)<=eps & xyzbDC(topologybDC(:,3),1)<=eps);
ElDCr=find( xyzbDC(topologybDC(:,1),1)>=lx-eps & xyzbDC(topologybDC(:,2),1)>=lx-eps & xyzbDC(topologybDC(:,3),1)>=lx-eps);
BCelemDC=ones(size(topologybDC,1),1);
BCelemDC(ElDCe,1)=2; BCelemDC(ElDCr,1)=3;
[BCtopo_DC,BCnodeA_DC,BCnodeB_DC]=bound(xyzbDC,topologybDC,BCelemDC,'n'); % Creates translation matrices for the node splitting accounting
MdcS=size(BCnodeB_DC,1);        


% Plot geometry
plotsurf(xyzbDC,topologybDC,[ElDCe;ElDCr],[NodDCe;NodDCr]); title(['Rect. duct - Nodes: ' num2str(size(xyzbDC,1)) ' - Elements: ' num2str(size(topologybDC,1))]); axis equal


% Prescribed normal velocity vector
if SplitNodes
    vnDC=zeros(size(BCnodeB_DC,1),1);
    vnDC(BCnodeB_DC(:,2)==2,1)=va;
else
    vnDC=zeros(Mdc,1);
    vnDC(NodDCe,1)=va;
end

% Impedance definition
if SplitNodes
    Ytmp=zeros(size(BCnodeB_DC,1),1);
    if ImpE
        Ytmp(BCnodeB_DC(:,2)==2 | BCnodeB_DC(:,2)==3,1)=EndAdm; % rho*c impedance is set at both ends, emitting and receiving
        vnDC=vnDC*2; % An anechoic piston generates half the pressure. Here it is corrected to have 1 Pa.
    else
        Ytmp(BCnodeB_DC(:,2)==3,1)=EndAdm; % rho*c impedance is set receiving end
    end
    YvectS=zeros(size(BCnodeB_DC,1),Mdc);
    for yy=1:length(Ytmp) % Create row-expanded admittance matrix
        YvectS(yy,BCnodeB_DC(yy,1))=Ytmp(yy);
    end
else
    Ytmp=zeros(Mdc,1);
    if ImpE
        Ytmp([NodDCe;NodDCr],1)=EndAdm; % rho*c impedance is set at both ends, emitting and receiving
        vnDC=vnDC*2; % An anechoic piston generates half the pressure. Here it is corrected to have 1 Pa.
    else
        Ytmp(NodDCr,1)=EndAdm; % rho*c impedance is set receiving end
    end
    YvectS=diag(Ytmp);
end
clear Ytmp


% Calculations
% Rectangular duct
paDC=zeros(Mdc,length(f)); % The acoustic pressure for each frequency
ptDC=zeros(length(f),1);  % Averaged pressure on the receiving end

for ii=1:length(f)
    disp(['Processing f= ' num2str(f(ii)) ' Hz'])
    
    if exist([fileGEOMduct '_split_f' num2str(floor(f(ii))) '.mat'],'file')==2
        load([fileGEOMduct '_split_f' num2str(floor(f(ii)))],'Ad','Bd','Cd');
    else
        [Ad,Bd,Cd]=TriQuadEquat(xyzbDC,[topologybDC(:,1:end-1) topologybDC(:,end)+j*BCelemDC],kp(ii)); % 3DBEM Version 4
        save([fileGEOMduct '_split_f' num2str(floor(f(ii)))],'Ad','Bd','Cd');
    end
    
    if ~SplitNodes
        BdS=Bd; Bd=zeros(Mdc);for bb=1:MdcS, Bd(:,BCnodeB_DC(bb,1))=Bd(:,BCnodeB_DC(bb,1))+BdS(:,bb);end; clear BdS % obtain unsplit version from the split Bd
    end
    
    % Solve system
    paDC(:,ii) = (Ad+j*kp(ii)*rho*c*Bd*YvectS) \ (-j*kp(ii)*rho*c*Bd*vnDC);
    ptDC(ii) = mean(abs(paDC(NodDCr,ii)));
    clear Ad Bd
    
end
clear YvectS


% Show pressure amplitude over length
figure;
jj=1; % Frequency to be plotted
plot(xyzbDC(:,1),abs(paDC(:,jj)),'x');grid
axis([0 lx 0 1.1]);
xlabel('Duct length');ylabel('|p| (Pa)');title(['Freq.: ' num2str(f(jj)) ' Hz, Std. wave ratio, p(max)/p(min): ' num2str(max(abs(paDC(:,1)))/min(abs(paDC(:,1))))]);

 
% Animated surface pressure
jj=1; % Frequency to be plotted
result=paDC(:,jj);
ap=[-30 40]; cca=[-max(abs(result)) max(abs(result))];
figure;  % screen animation
phase=0;
while 1 % Execution may be stopped using Ctr-C
    clf;%figure(hh)
    phase=phase+pi/180*10;    if phase>=2*pi-eps, phase=0;end;
    patch('faces',topologybDC(:,1:(end-1)),'vertices',xyzbDC(:,1:(end-1)),'FaceVertexCData',real(result*exp(j*phase)),'FaceColor','interp');
    title(['Rect. duct, ' num2str(f(jj)) ' Hz, end p. ' num2str(ptDC(jj,:)) ' Pa']);
    colorbar vert; caxis(cca);view(ap);axis equal
    drawnow;  %   pause(1)
end
close(hh)

 

% Animated pressure over length
jj=1; % Frequency to be plotted
figure;  % screen animation
phase=0;
result=paDC(:,jj);
while 1 % Execution may be stopped using Ctr-C
    clf;%figure(hh)
    phase=phase+pi/180*10;    if phase>=2*pi-eps, phase=0;end;
    plot(xyzbDC(:,1)*ones(1,length(f)),real(result*exp(j*phase)),'x');grid
    axis([min(xyzbDC(:,1)) max(xyzbDC(:,1)) -1.1 1.1]);
    title(['Rect. duct, ' num2str(f(jj)) ' Hz, end p. ' num2str(ptDC(jj,:)) ' Pa, lossy']);
    drawnow;  %   pause(1)
end
close(hh)
