% AxiBEM example: Calculates the scattering of plane wave by a sphere
% Compare with "FEM_ScatSphere_Axi_NuAc.m" in axisymmetrical FEM, Matlab PDE toolbox, (Numerical Acoustics) 

clear

R=1;               % Radius of the sphere
m=0;               % Excitación ejesimétrica, modo circunferencial m=0
quadelem=1;        % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)


% AMBIENT CONDITIONS
pa = 101325;          % Atmosferic pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 


% Wavenumber and frequency
% fr=150; 
% k=2*pi*fr/c0;
k=5;
fr=k*c/(2*pi); 
el_wl=6*max(fr)/c;  % Minimum mesh density

 

% GENERATION OF THE DOMAIN GEOMETRY
segments=[0 R R 0 20 R 0; R 0 0 -R 20 R 0];
%segments=[0 R 0 -R 20 R 0];
[rzb,topology]=nodegen(segments,'y',{},quadelem); % nodes and elements
M=size(rzb,1);N=size(topology,1);     % M nodos, N elementos

% CHIEF points:
rzb_chief=[linspace(0.1*R,0.55*R,10)' linspace(-0.8*R,0.6*R,10)' ones(10,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md')

% <<<<<<<<<<<<< Run to this point to see the geometry


% BEM MATRICES CALCULATION

% Incident sound field
%pI=[exp(i*k*rzb(:,2)) ; exp(i*k*rzb_chief(:,2))];
sources=[0 0 1 0]; % using the "incoming" function
pI=incoming(sources,[rzb ; rzb_chief],k,m);

% BEM calculation
[A,B,C]=BEMEquat0(rzb,topology,k,m,rzb_chief);

% Solution calculation on the boundary
ps=A\(-4*pi*pI);
vp=zeros(M,1); % Velocity is equal to 0 along the cylinder surfaces


% FIELD POINTS

% Field point definition and calculation using a mesh
espac=1/el_wl;
sizeFP=2*R;
nFP=max(sizeFP*el_wl,30);
rr=linspace(0,sizeFP,nFP); Nrr=length(rr);
zz=linspace(-sizeFP,sizeFP,nFP); Nzz=length(zz); 
[RR,ZZ]=meshgrid(rr,zz); % Mesh definition on the rho-z plane

[p_field,RRb,ZZb]=FieldPnt(RR,ZZ,ps,vp,rzb,topology,k,m,rho,c); % See FieldPnt help text

% Incident sound field on the field points:
%pIfp=exp(i*k*ZZb);
pIfp=incoming(sources,RRb,ZZb,k,m);


% POSTPROCESSING


% Analytical solution on the boundary
pAN = PlaneWaveScatSphere(k,R,R,linspace(0,pi,M),1e-6); pAN=conj(pAN);


% plot the pressure on the surface
figure;
thetaS = angle(rzb(:,1) + j*rzb(:,2))*180/pi;
subplot(2,1,1);
plot(thetaS,abs(ps),'ko--',thetaS,abs(pAN).','kx-');
title(['Scattering by a sphere - Frequency = ' num2str(fr) ' Hz']);
xlabel('Nodes on the generator'); ylabel('Pressure modulus (dB)')
axx=axis; axis([min(thetaS) max(thetaS) axx(3:4)])
legend('BEM','Analytical');grid;
subplot(2,1,2);
plot(thetaS,angle(ps)*180/pi,'ko--',thetaS,angle(pAN).'*180/pi,'kx-');
xlabel('Nodes on the generator'); ylabel('Phase of the pressure (degrees)')
axx=axis; axis([min(thetaS) max(thetaS) axx(3:4)])
legend('BEM','Analytical');grid;



% Analytical solution FP:
rFP = abs(RRb(1:end) + j*ZZb(1:end));
thetaFP = angle(RRb(1:end) + j*ZZb(1:end));
pfFPAN=zeros(length(rFP),1);
for fp=1:length(rFP)
    pfFPAN(fp,1) = PlaneWaveScatSphere(k,R,rFP(fp),thetaFP(fp)-pi/2,1e-6);
end
pfFPAN = conj(pfFPAN);


% plot the result on the field points, BEM and analytical
[r_i,c_i]=find(abs(RRb+j*ZZb)<=R); % remove results inside the cylinder
i_int=find(abs(RRb(1:end)+j*ZZb(1:end))<=R); 
pfFPAN(i_int)=NaN;
p_field_TOT = p_field + pIfp;
for ll=1:length(r_i)
    p_field_TOT(r_i(ll),c_i(ll))=NaN; 
end
figure;
subplot(2,1,1)
surf(RRb,ZZb,abs(p_field_TOT));
xlabel('r, m'); ylabel('z, m');zlabel('|p|, Pa');
title(['Scattering by a sphere (BEM), ka = ' num2str(k*R)])
axis equal
subplot(2,1,2)
surf(RRb,ZZb,reshape(abs(pfFPAN),size(RRb)));
xlabel('r, m'); ylabel('z, m');zlabel('|p|, Pa');
title(['Scattering by a sphere (Analytical), ka = ' num2str(k*R)])
axis equal
rotate3d on;



% ANIMATED results, displacement as z value
figure
hh=figure;  % screen animation
ph=0; AziEle=[10,80]; MaxMin=[-2 +2];
while 1 % Execution may be stopped by closing the figure or using Ctr-C
    figure(hh)
    ph=ph+pi/180*10;    if ph>=2*pi-eps, ph=0;end; % sweep over angles
    subplot(1,2,1)
    surf(RR,ZZ,real(p_field_TOT*exp(j*ph)));
    xlabel('r, m'); ylabel('z, m');zlabel('Re[p_{total}], Pa');
    title(['Sphere (BEM), ka = ' num2str(k*R)])
    axis equal
    axx=axis; axis([axx(1:4) MaxMin])
    rotate3d on; view(AziEle)
    subplot(1,2,2)
    surf(RR,ZZ,reshape(real(pfFPAN*exp(j*ph)),size(RR)));
    xlabel('r, m'); ylabel('z, m');zlabel('Re[p_{total}], Pa');
    title(['Sphere (Analytical), ka = ' num2str(k*R)])
    axis equal
    axx=axis; axis([axx(1:4) MaxMin])
    rotate3d on; view(AziEle)
    drawnow;  %   pause(1)
end


