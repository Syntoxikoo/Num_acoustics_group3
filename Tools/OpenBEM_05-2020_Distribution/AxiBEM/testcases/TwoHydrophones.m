% Two hydrophones face to face with plane waves coming from an angle
clear

% Geometrical parameters:
Rh=0.1;                 % Radius of the hydrophone, m
Lh=0.5;                 % Length of the hydrophone, m
Dh=0.1;                 % Distance face to face, m
A0=1;                   % Wave amplitude

% CONSTANTS
rho= 1026;              % kg/m3 density, sea water 13 C (Kinsler)
c= 1500;                % m/s speed of sound, sea water 13 C (Kinsler)

% General parameters:
fr=1000;               % Frequency (Hz)
kp=2*pi*fr/c;          % Wavenumber (1/m)
el_wl=6*max(fr)/c;     % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;         % Spacing of field points
AngleS=45*(pi/180);    % Angle of incidence of a plane wave in the rho-z plane, radians
phi_angle=0*(pi/180);  % Angle for the visualization plane in the circumferential direction, radians 
mmax=5;                % Number of terms in the cosine expansion (non-axisymmetrical conditions, see PhD thesis by Peter M Juhl)


% GENERATION OF THE DOMAIN GEOMETRY
segments=[0 Dh/2+Lh Rh Dh/2+Lh 3 0 el_wl;
          Rh Dh/2+Lh Rh Dh/2 15 0 el_wl;
          Rh Dh/2 0 Dh/2 3 0 el_wl;
          0 -Dh/2 Rh -Dh/2 3 0 el_wl;
          Rh -Dh/2 Rh -Dh/2-Lh 15 0 el_wl;
          Rh -Dh/2-Lh 0 -Dh/2-Lh 3 0 el_wl];
[rzb,topology]=nodegen(segments,'y'); % nodes and elements
M=size(rzb,1);N=size(topology,1);     % M nodos, N elementos

% CHIEF points inside the interior domain:
rzb_chief=[linspace(Rh*0.1,Rh*0.8,5)' linspace(Dh/2+Lh*0.11,Dh/2+Lh*0.91,5)' ones(5,1);
           linspace(Rh*0.1,Rh*0.8,5)' linspace(-Dh/2-Lh*0.11,-Dh/2-Lh*0.91,5)' ones(5,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md')


% <<<<<<<<<<<<< Run to this point to see the geometry

psm=zeros(M,mmax);  % The cosine coefficents for every m term forming the pressure on the surface is stored in this matrix
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

% Figure: pressure modulus over the surface, on a plane at angle "phi_angle" from the rho-z plane
figure;
ps_fb=[ps_front; ps_back];
plot3([rzb(:,1); -rzb(:,1)],[rzb(:,2); rzb(:,2)],abs(ps_fb),'*');grid
title(['|p| (Pa), freq.= ' num2str(fr) ' Hz, Inc angle.= ' num2str(AngleS/pi*180) ' deg., Graph angle: ' num2str(phi_angle/pi*180) ' deg.']);xlabel('r (m)');ylabel('z (m)')
axis equal
rotate3d on
    

% FIELD POINTS %%%%%%%%%%%%%%%%%%%

% Field points definition and calculation using a mesh
Nrr=10; rr=linspace(0,2*Rh,Nrr);
Nzz=10; zz=linspace(-Dh/2*0.9,Dh/2*0.9,Nzz);%linspace(-Dh-Lh,Dh+Lh,Nzz);
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
subplot(2,1,1)
contourf([-RR(:,end:-1:1) RR],[ZZ ZZ],abs([p_field_back(:,end:-1:1) p_field_front]));
title(['|p| (Pa), freq.= ' num2str(fr) ' Hz, Inc angle.= ' num2str(AngleS/pi*180) ' deg., Graph angle: ' num2str(phi_angle/pi*180) ' deg.']);xlabel('r (m)');ylabel('z (m)')
colorbar
subplot(2,1,2)
contourf([-RR(:,end:-1:1) RR],[ZZ ZZ],angle([p_field_back(:,end:-1:1) p_field_front])*180/pi);
title('p phase (degrees)')
colorbar



