% Sphere oscillating along the z-axis with viscous losses. There is an
% analytical solution in "Elements of Acoustics" by Temkin, sec. 6.9, and
% implemented in the Matlab function "SphereFirstOrder.m"

% Version 2: Assembled version including field points and figures. Used in 31265 Numerical Acoustics

clear, tic

% Input:
u0=1e-2;        % Amplitude of the oscillation of the sphere
f=200;ii=1;     % Frequency
Rad=1;%0.01;    % Sphere radius
N=20;%64;       % Number of elements

m=0;

% Constants
[rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu]=VTconst(f,-1);

% Create mesh and physical parameters
segments=[0 Rad 0 -Rad N Rad 0];

[rzb,topology,rzline]=nodegen(segments,'n'); 
M=size(rzb,1);     % M nodes, N elements
N=size(topology,1); 

[nvect,tvect]=normals(rzb,topology,'y');

% CHIEF points:
rzb_chief=[linspace(0.05,0.6,7)'*Rad linspace(-0.7,0.5,7)'*Rad ones(7,1)];
hold on; plot(rzb_chief(:,1),rzb_chief(:,2),'md')

% Constant matrices used in the calculation, a function of normal and tangential vectors:
NT1=nvect(:,1)*nvect(:,1)'+nvect(:,2)*nvect(:,2)'; NT2=nvect(:,1)*tvect(:,1)'+nvect(:,2)*tvect(:,2)';
NT3=tvect(:,1)*nvect(:,1)'+tvect(:,2)*nvect(:,2)'; NT4=tvect(:,1)*tvect(:,1)'+tvect(:,2)*tvect(:,2)';

% Gradient and laplacian matrices including the whole generator and varying spacings. Based on the secant formula applied twice.
[D0,Dl] = genderiv3(rzb,rzline);
% [D0,Dl] = genderiv6(rzb,topology,rzline);

% Normal and tangential velocity (first-order oscillating sphere)
vn0=u0*rzb(:,2)/Rad;
vt0=-u0*rzb(:,1)/Rad;

CondNrs=[]; % Condition numbers

% Calculate coefficient matrices
[Aa,Ba,CConsta]=BEMEquat0(rzb,topology,ka(ii),m,rzb_chief);
[Ah,Bh,CConsth]=BEMEquat0(rzb,topology,kh(ii),m,rzb_chief);
[Av,Bv,CConstv]=BEMEquat0(rzb,topology,kv(ii),m,rzb_chief);
%CondNrs=[CondNrs; cond(Aa) cond(Ba) cond(Ah) cond(Bh) cond(Av) cond(Bv)];
CondNrs=[CondNrs; cond(Aa) cond(Ba) cond(Av) cond(Bv)];

% Calculations
pas=zeros(M,1);     % The acoustic pressure
phs=zeros(M,1);     % The thermal pressure
pss=zeros(M,1);     % The lossless pressure

% Matrix multiplications
C1=phi_a(ii) - phi_h(ii)*tau_a(ii)/tau_h(ii);

% Include thermal contribution
ABah=phi_a(ii)*(Ba\Aa) - phi_h(ii)*(Bh\Ah)*tau_a(ii)/tau_h(ii);

% Include viscothermal contribution
rho_i_0=find(abs(rzb(:,1))<1e-10); % the values where rho=0 are changed to some very small number to avoid 1/r singularity in the divergence
rho0=rzb(:,1);
rho0(rho_i_0)=1e-10; %%%%% This value should not be too small (instability of the system) nor too high (unrealistic geometry). Work it out as a function of the other rho values
    
% New version of the previous corrected vectorial relationship. This one has been checked for equivalence in a dummy test
ABv=C1*inv(NT1.*(Bv\Av) + diag(1./rho0).*diag(nvect(:,1))) * ( NT2.*(Bv\Av)*D0 + Dl + diag(1./rho0).*diag(tvect(:,1)).*D0); 

%C2=-inv((rho0*ones(1,M)).*NT1.*(Bv\Av) + diag(nvect(:,1)))*((rho0*ones(1,M)).*(NT2.*(Bv\Av)) + ((rho0*ones(1,M)).*D0) + diag(tvect(:,1)));
C2=   -inv(NT1.*(Bv\Av) + diag(1./rho0).*diag(nvect(:,1))) * ( NT2.*(Bv\Av) + D0 +  diag(1./rho0).*(diag(tvect(:,1))));

rhs=vn0 - C2*vt0;
lhs=ABah + ABv;

% Solve system of equations:
pas=lhs\rhs;

% Local components of the viscous velocity:
v_r=ABv*pas + C2*vt0;
v_theta=vt0-(C1*D0)*pas;


% Analytical solution on the boundary and comparison plot
[pasAN, v_rAN, v_thetaAN, v_rAN_A, v_thetaAN_A, v_rAN_V, v_thetaAN_V] = SphereFirstOrder(ka,c,rho,Rad,u0,[Rad*ones(M,1) acos(rzb(:,2)/Rad)],-1,kv);
figure;
subplot(3,1,1)
plot(acos(rzb(:,2)/Rad)*180/pi,abs(pasAN),'-ko',acos(rzb(:,2)/Rad)*180/pi,abs(pas),'-rx'); grid
xlabel('Angle');ylabel('|pressure|');
title(['Acoustic variables on the surface, f = ' num2str(f(ii)) ' Hz, ' num2str(M) ' Nodes, sphere radius ' num2str(Rad)]);
legend('Analytical','BEM')
subplot(3,1,2)
plot(acos(rzb(:,2)/Rad)*180/pi,abs(v_rAN_V),'-ko',acos(rzb(:,2)/Rad)*180/pi,abs(v_r),'-rx'); grid
xlabel('Angle');ylabel('|Normal viscous velocity|');
legend('Analytical','BEM')
subplot(3,1,3)
plot(acos(rzb(:,2)/Rad)*180/pi,abs(v_thetaAN_A),'-ko',acos(rzb(:,2)/Rad)*180/pi,abs((C1*D0)*pas),'-rx'); grid
xlabel('Angle');ylabel('|Tangential acoustic velocity|');
legend('Analytical','BEM')



%% FIELD POINTS - FIELD POINTS - FIELD POINTS - FIELD POINTS - FIELD POINTS - FIELD POINTS 

% Acoustic velocity on the boundary
phs=-pas*tau_a(ii)/tau_h(ii); 
va=(Ba\Aa)*pas;
vh=(Bh\Ah)*phs; 

% Viscous velocity components and normal derivatives on the boundary in local coordinates
vt=v_theta;
vn=v_r; d
dvndn=NT1.*(Bv\Av)*vn + NT2.*(Bv\Av)*vt;
dvtdn=NT3.*(Bv\Av)*vn + NT4.*(Bv\Av)*vt;

% Viscous velocity components and normal derivatives on the boundary in global coordinates
vrho=nvect(:,1).*vn + tvect(:,1).*vt;
vz  =nvect(:,2).*vn + tvect(:,2).*vt;
dvrhodn=nvect(:,1).*dvndn + tvect(:,1).*dvtdn;
dvzdn  =nvect(:,2).*dvndn + tvect(:,2).*dvtdn;


% FIELD POINTS DEFINITION AND CALCULATION
dv=2.1e-3/sqrt(f); % Boundary layer thickness
% Grid 1: Line normal to the equator
Nrr1=20;rr1=linspace(Rad+1e4*eps,Rad+2*dv,Nrr1);

% Create two very close lines of field points in order to take gradient
ddv=dv/Nrr1/20; % separation distance between points in the direction perpendicular to the radius.

% The angle of the line of field points can be chosen
thetaFP=pi/2;%4; 
RR1=[rr1*sin(thetaFP) + ddv/2*cos(thetaFP); rr1*sin(thetaFP); rr1*sin(thetaFP) - ddv/2*cos(thetaFP)];
ZZ1=[rr1*cos(thetaFP) + ddv/2*sin(thetaFP); rr1*cos(thetaFP); rr1*cos(thetaFP) - ddv/2*sin(thetaFP)];

[vrho_field1,RRb1,ZZb1,Cfield_vr]=FieldPnt(RR1,ZZ1,vrho,dvrhodn,rzb,topology,kv(ii),m,0,c,0);
[vz_field1,RRb1,ZZb1,Cfield_vz]=FieldPnt(RR1,ZZ1,vz,dvzdn,rzb,topology,kv(ii),m,0,c,0);
[pa_field1,RRb1,ZZb1,Cfield_pa]=FieldPnt(RR1,ZZ1,pas,va,rzb,topology,ka(ii),m,0,c,0);
[ph_field1,RRb1,ZZb1]=FieldPnt(RR1,ZZ1,phs,vh,rzb,topology,kh(ii),m,0,c,0);
p_field1=pa_field1+ph_field1;

vr_field1     = vrho_field1*sin(thetaFP) + vz_field1*cos(thetaFP);
vtheta_field1 = vrho_field1*cos(thetaFP) - vz_field1*sin(thetaFP);


% PARTICLE VELOCITY CALCULATION
espacR=diff(rr1(1:2));espacZ=ddv;
[GRar,GRaz]=gradient(pa_field1([1 3],:),espacR,espacZ); 
[GRhr,GRhz]=gradient(ph_field1([1 3],:),espacR,espacZ); 
vr1=phi_a(ii)*GRar+phi_h(ii)*GRhr;
vz1=phi_a(ii)*GRaz+phi_h(ii)*GRhz;


%% DRAW FIGURES


% VISCOUS tangential velocity (several phases)
% To be used in PAPER. This figure is equivalent to figure 10-2 in Pierce, p. 526:
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
figure;
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
phase=[0:45:350]*pi/180;
%plot(raxis,abs(v_thetaAN_V_field1),'-ko',raxis,abs(vtheta_field1),'-rx'); grid
plot(raxis,real(-v_thetaAN_V_field1*exp(j*phase(1))),'-k',raxis,real(vtheta_field1(2,:)*exp(j*phase(1))),'r*'); grid; hold on
for pp=2:length(phase)
    plot(raxis,real(-v_thetaAN_V_field1*exp(j*phase(pp))),'-k',raxis,real(vtheta_field1(2,:)*exp(j*phase(pp))),'r*');
end
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential viscous velocity (m/s)','fontsize', 18);
%axis([0 2.1e-3/sqrt(f) -1e-2 1e-2])
axis([0 max(raxis) -u0*1.1 u0*1.1])
%legend('Analytical','BEM')
set(gca,'fontsize', 18)
%print -deps2 Temkin1


% TOTAL tangential velocity (several phases)
% To be used in PAPER. This figure is equivalent to figure 10-2 in Pierce, p. 526:
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
figure;
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
phase=[0:45:350]*pi/180;
%plot(raxis,abs(v_thetaAN_V_field1),'-ko',raxis,abs(vtheta_field1),'-rx'); grid
plot(raxis,real(-v_thetaAN_field1*exp(j*phase(1))),'-k',raxis,real((vz1(1,:)+vtheta_field1(2,:))*exp(j*phase(1))),'r*'); grid; hold on
for pp=2:length(phase)
    plot(raxis,real(-v_thetaAN_field1*exp(j*phase(pp))),'-k',raxis,real((vz1(1,:)+vtheta_field1(2,:))*exp(j*phase(pp))),'r*');
end
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential velocity (m/s)','fontsize', 18);
axis([0 max(raxis) -u0*1.1 u0*1.1])
set(gca,'fontsize', 18)
%print -deps2 Temkin1


% VISCOUS tangential velocity (modulus and phase)
% To be used in PAPER. This figure is equivalent to figure 10-2 in Pierce, p. 526:
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
figure;
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
subplot(2,1,1)
plot(raxis,abs(-v_thetaAN_V_field1),'-k',raxis,abs(vtheta_field1(2,:)),'r*'); grid; hold on
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential viscous velocity, modulus (m/s)','fontsize', 18);
axis([0 max(raxis) 0 u0*1.2])
set(gca,'fontsize', 18)
subplot(2,1,2)
plot(raxis,unwrap(angle(-v_thetaAN_V_field1))*180/2,'-k',raxis,unwrap(angle(vtheta_field1(2,:)))*180/2,'r*'); grid; hold on
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential viscous velocity, phase (deg.)','fontsize', 18);
axis([0 max(raxis) 0 360])
set(gca,'fontsize', 18)
%print -deps2 Temkin1


% ACOUSTIC+THERMAL tangential velocity (modulus and phase)
% To be used in PAPER. This figure is equivalent to figure 10-2 in Pierce, p. 526:
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
figure;
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
subplot(2,1,1)
plot(raxis,abs(-v_thetaAN_A_field1),'-k',raxis,abs(vz1(1,:)),'r*'); grid; hold on
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential acoustic velocity, modulus (m/s)','fontsize', 18);
axis([0 max(raxis) 0 u0*1.001])
set(gca,'fontsize', 18)
subplot(2,1,2)
plot(raxis,unwrap(angle(-v_thetaAN_A_field1))*180/2,'-k',raxis,unwrap(angle(-vz1(1,:)))*180/2,'r*'); grid; hold on
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential acoustic velocity, phase (deg.)','fontsize', 18);
axis([0 max(raxis) -360 1])
set(gca,'fontsize', 18)
%print -deps2 Temkin1


% TOTAL tangential velocity (modulus and phase)
% To be used in PAPER. This figure is equivalent to figure 10-2 in Pierce, p. 526:
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
figure;
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
subplot(2,1,1)
plot(raxis,abs(-v_thetaAN_field1),'-k',raxis,abs(vz1(1,:)+vtheta_field1(2,:)),'r*'); grid; hold on
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential velocity, modulus (m/s)','fontsize', 18);
axis([0 max(raxis) 0 u0*1.001])
set(gca,'fontsize', 18)
subplot(2,1,2)
plot(raxis,unwrap(angle(-v_thetaAN_field1))*180/2,'-k',raxis,unwrap(angle(-vz1(1,:)-vtheta_field1(2,:)))*180/2,'r*'); grid; hold on
xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential velocity, phase (deg.)','fontsize', 18);
axis([0 max(raxis) -360 1])
set(gca,'fontsize', 18)
%print -deps2 Temkin1

%% Animated plots

% VISCOUS tangential velocity (screen animated plot)
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
figure;hh=figure  % screen animation
phase=0;
while 1 % Execution may be stopped using Ctr-C
    clf;%figure(hh)
    phase=phase+pi/180*10;    if phase>=2*pi-eps, phase=0;end;
    plot(raxis,real(-v_thetaAN_V_field1*exp(j*phase)),'-k',raxis,real(vtheta_field1(2,:)*exp(j*phase)),'r*'); grid; hold on
    
    xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential viscous velocity (m/s)','fontsize', 18);
    axis([0 max(raxis) -u0*1.1 u0*1.1])
    set(gca,'fontsize', 18)
    title("VISCOUS tangential velocity")
    
    
    drawnow;  %   pause(1)
end
close(hh)
%%

% ACOUSTIC+THERMAL tangential velocity (screen animated plot)
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
figure;hh=figure  % screen animation
phase=0;
while 1 % Execution may be stopped using Ctr-C
    clf;%figure(hh)
    phase=phase+pi/180*10;    if phase>=2*pi-eps, phase=0;end;
    plot(raxis,real(-v_thetaAN_A_field1*exp(j*phase)),'-k',raxis,real((vz1(1,:))*exp(j*phase)),'r*'); grid; hold on
    
    xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential acoustic velocity (m/s)','fontsize', 18);
    axis([0 max(raxis) -u0*1.1 u0*1.1])
    set(gca,'fontsize', 18)
    title("ACOUSTIC+THERMAL tangential velocity")
    
    drawnow;  %   pause(1)
end
close(hh)

%%

% TOTAL tangential velocity (screen animated plot)
[pasAN_field1, v_rAN_field1, v_thetaAN_field1, v_rAN_A_field1, v_thetaAN_A_field1, v_rAN_V_field1, v_thetaAN_V_field1] = ...
    SphereFirstOrder(ka,c,rho,Rad,u0,[sqrt(RR1(2,:).^2+ZZ1(2,:).^2).' ones(size(RR1(2,:)))'*(thetaFP)],-1,kv);
raxis=sqrt(RR1(2,:).^2+ZZ1(2,:).^2)-Rad;
figure;hh=figure  % screen animation
phase=0;
while 1 % Execution may be stopped using Ctr-C
    clf;%figure(hh)
    phase=phase+pi/180*10;    if phase>=2*pi-eps, phase=0;end;
    plot(raxis,real(v_thetaAN_field1*exp(j*phase)),'-k',raxis,real((vz1(1,:)+vtheta_field1(2,:))*exp(j*phase)),'r*'); grid; hold on
    
    xlabel('Distance from the surface (m)','fontsize', 18);ylabel('Tangential velocity (m/s)','fontsize', 18);
    axis([0 max(raxis) -u0*1.1 u0*1.1])
    set(gca,'fontsize', 18)
    title("TOTAL tangential velocity")
    
    drawnow;  %   pause(1)
end
close(hh)





