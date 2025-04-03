clear

% ---------------------- PREPROCESSING ------------------------

% Ambient conditions
pa = 101325;         
t = 20;              
Hr = 50;             
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% Frequency
fr=1000; ff=1;
k=2*pi*fr/c;

% Source position
position = [0,1];

% nb elem
el_wl=6*max(fr)/c;  
espac=1/el_wl; 


% % Corresponding normal plane wave absorption and reflection index:
% alfa=4*real(1./betaBs)./(((abs(1./betaBs)).^2)/rho/c + 2*real(1./betaBs) + rho*c);
% R=((1./betaBs)-rho*c)./((1./betaBs)+rho*c);

N = 11;
lamb = 2*pi / k;
depths = lamb ./(2*N) *mod((0:N-1).^2,N);
point = (1:2*length(depths));
point(1:2:end) = depths;
point(2:2:end) = depths;


ext_width = 0.05;
int_width = 0.1;
pt1x = -N/2 * int_width - (N/2+1) * ext_width + (1:N)*(ext_width + int_width);
pt4x = -N/2 * int_width - (N/2) * ext_width + (1:N)*(ext_width + int_width);


% Create arrays for plotting vertical lines (wells)
x = [];
y = [];

% For each well
for i = 1:N
    % Add points to draw the well (in order)
    x = [x, pt1x(i), pt1x(i), pt4x(i), pt4x(i)];
    y = [y, 0, depths(i), depths(i), 0];

end
y = y+0.01;
y(1) = 0;y(end) = 0;
x = [x x(end) x(1)];
y = [y y(end) y(1)];

% Plot the diffuser shape
figure;
plot(x, y, 'b-', 'LineWidth', 1.5);
grid on;
title('Shape of Diffuser');
xlabel('X-axis');
ylabel('Y-axis');


% Define points for segments
point_1 = [x(1:end-1); y(1:end-1)]';
point_2 = [x(2:end); y(2:end)]';

% Calculate number of elements per segment based on element length
elem_num = max(ceil(sqrt(sum((point_2 - point_1).^2, 2)) * el_wl), 1);

% Create segments array
segmentsB = [point_1 point_2 elem_num zeros(size(point_1, 1), 2)];

% Generate nodes and topology using nodegen
[xybB, topologyB, rzline, segrzb] = nodegen(segmentsB, 'y');
[xybndB,topologyndB,xydumB,topodumB,xynodumB]=nodummy(xybB,topologyB,'y'); 

% Chieff point 

xyb_chief=[min(x) + (max(x)-min(x))*rand(10,1) 0.01 * rand(10,1) -ones(10,1)];

hold on;plot(xyb_chief(:,1),xyb_chief(:,2),'r^')


% % admittance of the barrier by segments
sigmaB=Inf(size(segmentsB,1),1); % Reflecting barrier  <<<<< UPDATE IF MORE SEGMENTS ARE INTRODUCED IN THE GEOMETRY 

betaBs=betag(sigmaB,fr)/(rho*c); % de-normalized admittance, as needed to solve the system and domain points

% plane adm 
sigmaP=Inf;
betaPs=betag(sigmaP,fr); 
admittB=betaBs(floor(segrzb+1/2));

[XX,YY]=meshgrid(-5:espac:5,0.01:espac:5);
xy=[XX(1:end)' YY(1:end)'];

% ---------------------- Calculation ------------------------
% obtain incident pressure on the mesh nodes and on the CHIEF points
[G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k(ff),position,betaPs(ff),... %xybndB(:,1),xybndB(:,2));
    [xybndB(:,1); xyb_chief(:,1)],[xybndB(:,2); xyb_chief(:,2)]);
inc_pressure=i/4*(G0dir+G0ref)+Pbeta;

% calculate coefficient matrices with CHIEF points
[Ab,Bb]=bem2d(xybB,topologyB,k(ff),betaPs(ff),xyb_chief);

% solve system
RHS=[-2*pi*inc_pressure]; 
pB=(Ab+1j*k*rho*c*Bb*diag(admittB(xynodumB)))\RHS;


[G0dir,G0ref,dG0dirdR1,dG0refdR2,Pbeta,dPbetadx,dPbetady]=greendef(k(ff),position,betaPs(ff),xy(:,1),xy(:,2));
pIfield=i/4*(G0dir+G0ref)+Pbeta;


% calculation of the pressure on field points:
[Ap,Bp,CConst]=fieldpoints(xybB,topologyB,k(ff),betaPs(ff),xy);
% solve the pressure on the field points and get the pressure relative to free field
pfield=((Ap+j*k*rho*c*Bp*diag(admittB(xynodumB)))*pB./CConst+pIfield).';

SPLff=20*log10(abs(pfield./pIfield.'));


% POSTPROCESSING

% plot the result on the field points
figure;
subplot(2,1,1);

plotfpoints(XX,YY,SPL(pfield)); hold on
% plot(xybndB(:,1), xybndB(:,2),LineWidth=2)
xlabel('Distance (m)');
ylabel('Height (m)');
title(['SPL (dB) - Frequency = ' num2str(fr(ff)) ' Hz']);
axis equal
subplot(2,1,2);
plotfpoints(XX,YY,SPLff); hold on
% plot(xybndB(:,1), xybndB(:,2),LineWidth=2)
xlabel('Distance (m)');
ylabel('Height (m)');
title(['SPL rel. free field (dB) - Frequency = ' num2str(fr(ff)) ' Hz']);
axis equal


