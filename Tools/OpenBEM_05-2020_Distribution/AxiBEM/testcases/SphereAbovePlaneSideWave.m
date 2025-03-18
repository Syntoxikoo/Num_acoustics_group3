% Sphere above a plane, with a plane wave coming from an angle
% Reflecting plane at z=0
% VCH, 5-2016

clear

% Geometrical parameters:
R=1;                    % Radius of the sphere
D=1.5*R;                % Distance between centre and plane
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
AngleS=90*(pi/180);    % Angle of incidence of a plane wave in the rho-z plane, radians (0 to 90, due to reflecting plane)
phi_angle=0*(pi/180);  % Angle for the visualization plane in the circumferential direction, radians 
mmax=8;                % Number of terms in the cosine expansion (non-axisymmetrical conditions, see PhD thesis by Peter M Juhl)
quadelem=1;            % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)
planeON=1;             % This example has a reflecting plane, this must be 1.

% GENERATION OF THE DOMAIN GEOMETRY
segments=[0 R+D R D 20 R el_wl; R D 0 -R+D 20 R el_wl];
[rzb,topology]=nodegen(segments,'y',{},quadelem); % nodes and elements
M=size(rzb,1);N=size(topology,1);     % M nodes, N elements


% CHIEF points inside the interior domain:
rzb_chief=[linspace(R*0.1,R*0.6,5)' linspace(R*0.11,R*0.5,5)'+D ones(5,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md')


% <<<<<<<<<<<<< Run to this point to see the geometry

% Obtain some of the variables needed for the modified on-plane incident
% pressure, solving the system and boundary conditions (velocity, admittance)
[rzbnd,topologynd,rzdum,topodum,rznodum]=nodummy(rzb,topology,'y'); 
Mnd=size(rzbnd,1);

psm=zeros(Mnd,mmax);  % The cosine coefficents for every m term forming the pressure on the surface is stored in this matrix
% Single plane wave at an angle. Several sources could be defined here, but
% the sum would be coherent. Note that the image source is also defined.
sources=[AngleS 0 A0 0 ; pi-AngleS 0 A0 0]; 

% BEM calculation of the surface pressure
for m=0:mmax-1
    % Incident sound field
    pI=incoming(sources,[rzbnd ; rzb_chief],kp,m);
    
    % BEM calculation
    A=BEMEquat0(rzb,topology,kp,m,rzb_chief,planeON);
    
    % Calculation on the boundary
    psm(:,m+1)=A\(-4*pi*pI);
end

% Perform the summation of the calculated m terms of the pressure for the given phi angle
ps_front=zeros(Mnd,1);
ps_back=zeros(Mnd,1);
for m=0:mmax-1
    ps_front=ps_front+psm(:,m+1)*cos(m*phi_angle);
    ps_back=ps_back+psm(:,m+1)*cos(m*(phi_angle+pi));
end
vp=zeros(Mnd,1);  % Velocity is equal to 0 along the boundary (hard surfaces)


% Figure: pressure modulus over the surface, on a plane at angle "phi_angle" from the rho-z plane
figure;
ps_fb=[ps_front;ps_back(end-1:-1:1)];
plot(linspace(0,360,2*Mnd-1),abs(ps_fb)/2,'-x');grid
title(['|p| (Pa), freq.= ' num2str(fr) ' Hz, Inc angle.= ' num2str(AngleS/pi*180) ' deg., Slice angle: ' num2str(phi_angle/pi*180) ' deg.']);
xlabel('theta (deg.)');ylabel('Amplitude modulus (Pa)')
legend('BEM','Location','southeast');


% FIELD POINTS %%%%%%%%%%%%%%%%%%%

% Field points definition and calculation using a mesh
Nrr=10; rr=linspace(0,2*R,Nrr);
Nzz=20; zz=linspace(0,5*R,Nzz);%linspace(-Dh-Lh,Dh+Lh,Nzz);
[RR,ZZ]=meshgrid(rr,zz); % Mesh definition on the rho-z plane

p_field_m=zeros(Nzz,Nrr,mmax);
for m=0:mmax-1
    % Incident sound field
    pI_field=incoming(sources,RR,ZZ,kp,m);
    
    p_scat_field=FieldPnt(RR,ZZ,psm(:,m+1),vp,rzb,topology,kp,m,rho,c,0,planeON); % See FieldPnt help text
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
% surf([-RR(:,end:-1:1) RR -RR(:,end:-1:1) RR],[ZZ ZZ -ZZ -ZZ],...
%     abs([p_field_back(:,end:-1:1) p_field_front p_field_back(:,end:-1:1) p_field_front])/2); % Represent mirror data and divide by 2
title(['|p| (Pa), freq.= ' num2str(fr) ' Hz, Inc angle.= ' num2str(AngleS/pi*180) ' deg., Slice angle: ' num2str(phi_angle/pi*180) ' deg.']);
xlabel('r (m)');ylabel('z (m)')
rotate3d on
colorbar



