%==================================================================
% Exercise 7.2
% Example: Impedance tube.
%
% Mesh: using quadratic elements
%
% Illustration of the procedure for the measurement of an impedance 
% of a sample using the two-microphone method
%
%==================================================================
clc; close all;
%% Input data
U0     = 1e-3;
rho0   = 1.2;   % Fluid density
c0     = 342.2; % Speed of sound
P0     = U0*rho0*c0^2;     % Input pressure
D      = 0.1;   % Piston diameter 
L      = 10*D;  % Length of the cavity
s      = D/2;   % microphone spacing
d      = D/2;   % distance between mic 2 and sample
%% Frequency domain
fc    = floor(1.84*c0/D/pi); % Cut off frequency
freq  = 100:2:fc;            % Choose correctly the lowest frequency
omega = 2*pi*freq;
k0    = omega/c0;
Nfreq = length(freq);
%% Sample properties
Z0    = rho0 * c0;
h     = 0.02;                % Thickness of the sample
sigma = 10000;               % Flow resitivity
X     = rho0*freq/sigma;
Zc    = Z0*(1+0.057*X.^(-0.754)-1i*0.087.*X.^(-0.732)); % characteristic impedance
k     = k0 .*(1+0.0978*X.^(-0.700)-1i*0.189.*X.^(-0.595));
Z     = -1i.*Zc.*cot(k*h) /Z0;
beta  = 1./Z;   % Convert to admittance
%% Step 1: Mesh
ne  = 30;                       % Number of quadratic elements (make sure that there are nodes at the two mic locations)
nnt = 2*ne+1;                   % Total number of nodes
h   = L/ne;                     % Length of the elements
x   = 0:h/2:L;                  % Coordinates table
in2 = find(abs(x-d)<1e-6,1);      % Location of mic 2
in1 = find(abs(x-(d+s))<1e-6,1);  % Location of mic 1
s   = abs(x(in1)-x(in2));       % recalculate microphones separation
d   = x(in2);                   % distance between mic 2 and the sample
%% Step 2: Compute Elementary matrices
Ke = [7,-8,1;-8,16,-8;1,-8,7]/(3*h)*c0^2;
Me = [4,2,-1;2,16,2;-1,2,4]*h/30;

%% Step 3: Assembling
I = eye(3,3);
K = zeros(nnt,nnt); M = zeros(nnt,nnt);

%% ------Write script here------
 % Task: assemble the globale matrices K and M using the quadratic elements Ke and Me.
 % Hint: use similar aproach as is found in 'FEM_ex7_1_eigenmodes.m', but using quadratic elements. 
 for ie = 1:ne
    idx = 2*ie-1;
    L = zeros(3,nnt); L(:,idx:idx+2) = I;
    K = K + L'* Ke*L;
    M = M + L'*Me*L;
 end

% K = K(2:nnt-1,2:nnt-1);
% M = M(2:nnt-1,2:nnt-1);
%% Step 4 & 5: specify frequency dependent impedance and Solve the system with the Force vector : Piston at a x=L
ndof         = nnt;
Pquad_direct = zeros(Nfreq,ndof);
Pquad_modal  = zeros(1,ndof);
Pquad_exact  = zeros(1,ndof);
P_mic1       = zeros(1,ndof);
P_mic2       = zeros(1,ndof);
A            = zeros(nnt,nnt);
f           = zeros(ndof,1); % Initialize input pressure vector

%% er
P_mic1 = zeros(1,Nfreq);
P_mic2 = zeros(1,Nfreq);

%% ------Write script here------
 % Loop for each frequency and solve FE equation S = K-w^2*M+i1*w*c0*A
 % P = S\f
 for ff = 1:Nfreq
    f(end,1) =P0;
    A(1,1) = beta(ii);
    
    S = K *omega(ff)^2 * M + 1j * omega(ff) * c0 *A;
    
    Pquad_direct(ff,:) = S \ f;
 end
P_mic1 = Pquad_direct(:,in1);
P_mic2 = Pquad_direct(:,in2);


% calculate the normalized impedance
k     = omega/c0;
H12   = P_mic2./P_mic1;
R     = (H12-exp(-1i*k*s))./(exp(1i*k*s)-H12) .*exp(1i*2*k*(d+s));
Z_num = (1+R)./(1-R);
%% Step 6: Comparison with the exact solution
figure(1)
subplot(1,2,1)
plot(freq,real(Z),'k','LineWidth',2)
hold on
plot(freq,real(Z_num),':','LineWidth',2,'Color',[0.5 0.5 0.5])
xlim([100 2000])
grid on
xlabel('Frequency (Hz)');
ylabel(' Normalized Impedance - Real part');
legend('Analytical', 'FEM ');
subplot(1,2,2)
plot(freq,imag(Z),'k','LineWidth',2)
hold on
plot(freq,imag(Z_num),':','LineWidth',2,'Color',[0.5 0.5 0.5]);
xlim([100 2000])
grid on
xlabel('Frequency (Hz)');
ylabel(' Normalized Impedance - Imaginary part');
figure (2)
R_theo=(Z-1)./(Z+1); alpha_theo=1-abs(R_theo).^2;
R_num=(Z_num-1)./(Z+1); alpha_num=1-abs(R_num).^2;
plot(freq,alpha_theo,'k','LineWidth',2)
hold on
plot(freq,alpha_num,':','LineWidth',2,'Color',[0.5 0.5 0.5])
xlim([100 2000])
grid on
xlabel('Frequency (Hz)');
ylabel(' Absorption coefficient');
legend('Analytical', 'FEM ');