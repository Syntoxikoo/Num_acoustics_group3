% Example of calculation: line source and barrier

clear

% PREPROCESSING

% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (�C)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% frequency and wavenumber
fr=150; ff=1; % (ff is the frequency count when a loop is used)
k=2*pi*fr/c;

% number of elements per wavelength
el_wl=6*max(fr)/c;   % Minimum mesh density as a function of the highest frequency
espac=1/el_wl;       % Spacing of field points

% admittance of the barrier by segments
sigmaB=[Inf Inf Inf Inf]; % Reflecting barrier  <<<<< UPDATE IF MORE SEGMENTS ARE INTRODUCED IN THE GEOMETRY 

betaBs=betag(sigmaB,fr)/(rho*c); % de-normalized admittance, as needed to solve the system and domain points

% betaBs=1/(rho*c)*ones(1,4); % set this, for example, to have a barrier impedance of rho*c

% Corresponding normal plane wave absorption and reflection index:
alfa=4*real(1./betaBs)./(((abs(1./betaBs)).^2)/rho/c + 2*real(1./betaBs) + rho*c);
R=((1./betaBs)-rho*c)./((1./betaBs)+rho*c);

% dimensions of the barrier
Bhgt=3; % barrier heigth
Bwdt=0.2;   % barrier width
bar=[-Bwdt 0;...
     -Bwdt Bhgt;...
      Bwdt/2 Bhgt;...
      Bwdt/2 0];
% barrier top
topsegs=[bar(2,:) bar(3,:) ceil(sqrt(sum((bar(2,:)-bar(3,:)).^2))*el_wl) 0 0]; % node 1 - node 2 - nbr of elements 

% topsegs=[bar(2,:) bar(3,:) 10 (0.5*Bwdt+eps) 0];   % Example another top

% TASK: Change the design of the barrier top (topsegs) and observe the result. Did
% you manage to improve the performance of the barrier?

% normalized plane admittance
sigmaP=Inf;
betaPs=betag(sigmaP,fr); % normalized admittance, as required by greendef, bem2d and fieldpoints

% line source position
position=[-10 0.5];

% generate mesh of field points
[XX,YY]=meshgrid(-11.25:espac:10,0.25:espac:10);
xy=[XX(1:end)' YY(1:end)'];


% define geometry
segmentsB=[bar(1,:) bar(2,:) ceil(sqrt(sum((bar(1,:)-bar(2,:)).^2))*el_wl) 0 0;... % Front side
    topsegs;...                                                                    % TOP
    bar(3,:) bar(4,:) ceil(sqrt(sum((bar(3,:)-bar(4,:)).^2))*el_wl) 0 0;...        % Back side
    bar(4,:) bar(1,:) ceil(sqrt(sum((bar(4,:)-bar(1,:)).^2))*el_wl) 0 0];          % Lower (ground) side

[xybB,topologyB,rzline,segrzb]=nodegen(segmentsB,'y');
admittB=betaBs(floor(segrzb+1/2));
[xybndB,topologyndB,xydumB,topodumB,xynodumB]=nodummy(xybB,topologyB,'y'); % say 'y' if there is a plane, otherwise 'n'.

% CHIEF points:
xyb_chief=[linspace(-Bwdt*0.2,Bwdt*0.3,10)' linspace(Bhgt*0.09,Bhgt*0.88,10)' -ones(10,1)];
hold on;plot(xyb_chief(:,1),xyb_chief(:,2),'md')

% CALCULATION
eps
% obtain incident pressure on the mesh nodes and on the CHIEF points
[G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k(ff),position,betaPs(ff),... %xybndB(:,1),xybndB(:,2));
    [xybndB(:,1); xyb_chief(:,1)],[xybndB(:,2); xyb_chief(:,2)]);
inc_pressure=i/4*(G0dir+G0ref)+Pbeta;

% calculate coefficient matrices with CHIEF points
[Ab,Bb]=bem2d(xybB,topologyB,k(ff),betaPs(ff),xyb_chief);

% solve system
RHS=[-2*pi*inc_pressure]; % Right-Hand-Side of the system
pB=(Ab+1j*k*rho*c*Bb*diag(admittB(xynodumB)))\RHS;

% obtain incident pressure on the field points
[G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k(ff),position,betaPs(ff),xy(:,1),xy(:,2));
pIfield=i/4*(G0dir+G0ref)+Pbeta;


% calculation of the pressure on field points:
[Ap,Bp,CConst]=fieldpoints(xybB,topologyB,k(ff),betaPs(ff),xy);
% solve the pressure on the field points and get the pressure relative to free field
pfield=((Ap+j*k*rho*c*Bp*diag(admittB(xynodumB)))*pB./CConst+pIfield).';

% % Field point calculation for cases with many field points (replaces the previous)
% for ii=1:size(xy,1);
%    [Ap,Bp,CConst]=fieldpoints(xybB,topologyB,k(ff),betaPs(ff),xy(ii,:));
%     pfield(ii)=((Ap+j*k*rho*c*Bp*diag(admittB(xynodumB)))*pB./CConst+pIfield(ii)).';
%    disp(['Field point nr. ' num2str(ii) ' of ' num2str(size(xy,1))])
%    clear Ap CConst
% end

SPLff=20*log10(abs(pfield./pIfield.'));


% POSTPROCESSING

% plot the result on the field points
figure;
subplot(2,1,1);

plotfpoints(XX,YY,SPL(pfield)); hold on
plot(xybndB(:,1), xybndB(:,2),LineWidth=2)
xlabel('Distance (m)');
ylabel('Height (m)');
title(['SPL (dB) - Frequency = ' num2str(fr(ff)) ' Hz']);
axis equal
subplot(2,1,2);
plotfpoints(XX,YY,SPLff); hold on
plot(xybndB(:,1), xybndB(:,2),LineWidth=2)
xlabel('Distance (m)');
ylabel('Height (m)');
title(['SPL rel. free field (dB) - Frequency = ' num2str(fr(ff)) ' Hz']);
axis equal

