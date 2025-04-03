% Sphere with a plane wave coming from an angle

% VCH, 5-2016
clear

% Geometrical parameters:
R=1;                    % Radius of the sphere
A0=1;                   % Wave amplitude

% AMBIENT CONDITIONS
pa = 101325;          % Atmosferic pressure (Pa)
t = 20;               % Temperature (ºC)
Hr = 50;              % Humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% General parameters:
fr=100;                % Frequency (Hz)
kp=2*pi*fr/c;          % Wavenumber (1/m)
el_wl=6*max(fr)/c;     % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;         % Spacing of field points
AngleS=0*(pi/180);    % Angle of incidence of a plane wave in the rho-z plane, radians
phi_angle=0*(pi/180);  % Angle for the visualization plane in the circumferential direction, radians 
mmax=8;                % Number of terms in the cosine expansion (non-axisymmetrical conditions, see PhD thesis by Peter M Juhl)
quadelem=0;            % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)


% GENERATION OF THE DOMAIN GEOMETRY
segments=[0 R R 0 20 R el_wl; R 0 0 -R 20 R el_wl];
[rzb,topology]=nodegen(segments,'y',{},quadelem); % nodes and elements
M=size(rzb,1);N=size(topology,1);     % M nodes, N elements


% CHIEF points inside the interior domain:
rzb_chief=[linspace(R*0.1,R*0.6,5)' linspace(R*0.11,-R*0.5,5)' ones(5,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md')


%% <<<<<<<<<<<<< Run to this point to see the geometry

psm=zeros(M,mmax);  % The cosine coefficents for every m term forming the pressure on the surface are stored in this matrix
sources=[AngleS 0 A0 0]; % Single plane wave at an angle. Several sources could be defined here, but the sum would be coherent.

% BEM calculation of the surface pressure
for m=0:mmax-1
    % Incident sound field
    pI=incoming(sources,[rzb ; rzb_chief],kp,m);
    
    % BEM calculation
    A=BEMEquat0(rzb,topology,kp,m,rzb_chief);
    
    % Calculation on the boundary
    psm(:,m+1)=A\(-4*pi*pI);
end

% Perform the summation of the calculated m terms of the pressure for the given phi angle
ps_front=zeros(M,1);
ps_back=zeros(M,1);
for m=0:mmax-1
    ps_front=ps_front+psm(:,m+1)*cos(m*phi_angle);
    ps_back=ps_back+psm(:,m+1)*cos(m*(phi_angle+pi));
end
vp=zeros(M,1);  % Velocity is equal to 0 along the boundary (hard surfaces)

% Analitical solution for verification:
ptotANA=PlaneWaveScatSphere(kp,R,R,linspace(-AngleS,-AngleS+2*pi,100),1e-6);

% Figure: pressure modulus over the surface, on a plane at angle "phi_angle" from the rho-z plane
figure;
ps_fb=[ps_front; ps_back(end-1:-1:1)];
plot(linspace(0,360,2*M-1),abs(ps_fb),'-x',linspace(0,360,length(ptotANA)),abs(ptotANA),'-');grid
title(['|p| (Pa), freq.= ' num2str(fr) ' Hz, Inc angle.= ' num2str(AngleS/pi*180) ' deg., Slice angle: ' num2str(phi_angle/pi*180) ' deg.']);
xlabel('theta (deg.)');ylabel('Amplitude modulus (Pa)')
legend('BEM','Analytical','Location','southwest');
    

% FIELD POINTS %%%%%%%%%%%%%%%%%%%

% Field points definition and calculation using a mesh
Nrr=10; rr=linspace(0,2*R,Nrr);
Nzz=20; zz=linspace(-2*R,2*R,Nzz);%linspace(-Dh-Lh,Dh+Lh,Nzz);
[RR,ZZ]=meshgrid(rr,zz); % Mesh definition on the rho-z plane

p_field_m=zeros(Nzz,Nrr,mmax);
for m=0:mmax-1
    % Incident sound field
    pI_field=incoming(sources,RR,ZZ,kp,m);
    
    p_scat_field=FieldPnt(RR,ZZ,psm(:,m+1),vp,rzb,topology,kp,m,rho,c,0); % See FieldPnt help text
    p_field_m(:,:,m+1)=p_scat_field + pI_field;
end

% Perform the summation of the calculated m terms of the field points pressure for a this phi angle
p_field_front=zeros(Nzz,Nrr);
p_field_back=zeros(Nzz,Nrr);
for m=0:mmax-1
    p_field_front=p_field_front + p_field_m(:,:,m+1)*cos(m*phi_angle);
    p_field_back=p_field_back + p_field_m(:,:,m+1)*cos(m*(phi_angle+pi));
end


% Plot sound pressure on the field point mesh
figure
surf([-RR(:,end:-1:1) RR],[ZZ ZZ],abs([p_field_back(:,end:-1:1) p_field_front]));
title(['|p| (Pa), freq.= ' num2str(fr) ' Hz, Inc angle.= ' num2str(AngleS/pi*180) ' deg., Slice angle: ' num2str(phi_angle/pi*180) ' deg.']);
xlabel('r (m)');ylabel('z (m)')
rotate3d on
colorbar



