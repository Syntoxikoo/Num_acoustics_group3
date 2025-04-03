
% Example of calculation: Infinite cylinder in free field
% CONVERGENCE TEST
clear

% PREPROCESSING

% Constants
rho=1.1992;
c=344;

% frequency and wavenumber
fr=150; 
k=2*pi*fr/c;



% Define numbers of elements per segment to be tried
NelT = round(logspace(1,2,10)/2);

CylError=ones(length(NelT),1);
for nn=1:length(NelT)
    
    Nel = NelT(nn); % Total number of elements /4
    
    % define geometry of the object
    Rc=1; % Radius of the cylinder
    segments=[-Rc 0 0 Rc Nel Rc 0;...
        0 Rc Rc 0 Nel Rc 0;...
        Rc 0 0 -Rc Nel Rc 0;...
        0 -Rc -Rc 0 Nel Rc 0];
    [xyb,topology]=nodegen(segments,'n');
    
    % CALCULATION
    
    % obtain incident pressure
    inc_pressure=exp(i*k*xyb(:,1));
    
    % calculate coefficient matrix
    A=bem2d(xyb,topology,k); 
    
    % solve system
    ps=A\(-2*pi*inc_pressure);
    
    % Analytical solution
    pAna=cylscat(k,Rc,xyb(:,1:2),150);
    
    CylError(nn) = mean(abs(ps-pAna.')./abs(pAna.'));
end


% POSTPROCESSING

% Convergence plot
figure;
loglog(NelT*4,CylError,'-ro');
title(['Convergence of BEM solution - k: ' num2str(k)])
xlabel('Number of elements'); ylabel('Relative error')
grid

%%

load("Ex_5/2DBEM_Exercise_5/k_fail.mat");

% Define numbers of elements per segment to be tried
NelT = round(logspace(1,2,10)/2);

CylError=ones(length(NelT),1);
for nn=1:length(NelT)
    
    Nel = NelT(nn); % Total number of elements /4

    
    % define geometry of the object
    Rc=1; % Radius of the cylinder
    segments=[-Rc 0 0 Rc Nel Rc 0;...
        0 Rc Rc 0 Nel Rc 0;...
        Rc 0 0 -Rc Nel Rc 0;...
        0 -Rc -Rc 0 Nel Rc 0];
    [xyb,topology]=nodegen(segments,'n');
    
    % CALCULATION
    
    % obtain incident pressure
    inc_pressure=exp(i*k*xyb(:,1));
    
    % calculate coefficient matrix
    A=bem2d(xyb,topology,k); 
    
    % solve system
    ps=A\(-2*pi*inc_pressure);
    
    % Analytical solution
    pAna=cylscat(k,Rc,xyb(:,1:2),150);
    
    CylError(nn) = mean(abs(ps-pAna.')./abs(pAna.'));
end

% POSTPROCESSING

% Convergence plot
figure;
loglog(NelT*4,CylError,'-ro');

title(['Convergence of BEM solution - k: ' num2str(k)])
xlabel('Number of elements'); ylabel('Relative error')
grid

%%
% Define numbers of elements per segment to be tried
NelT = round(logspace(1,2,10)/2);

ChiefPoint = [0.5,-0.2 ; -0.35,0.74];

xCh = -1 + 2* rand([4,1],"double");
yCh = -1+2*rand([4,1],"double");
ChiefPoint = [xCh,yCh];

CylError2=ones(length(NelT),1);
for nn=1:length(NelT)
    
    Nel = NelT(nn); % Total number of elements /4
    
    % define geometry of the object
    Rc=1; % Radius of the cylinder
    segments=[-Rc 0 0 Rc Nel Rc 0;...
        0 Rc Rc 0 Nel Rc 0;...
        Rc 0 0 -Rc Nel Rc 0;...
        0 -Rc -Rc 0 Nel Rc 0];
    [xyb,topology]=nodegen(segments,'n');
    
    % CALCULATION
    
    % obtain incident pressure
    inc_pressure=exp(i*k*xyb(:,1));
    pIch=exp(1j*k*ChiefPoint(:,1));
    
    % calculate coefficient matrix
    A=bem2d(xyb,topology,k); 
    [Apchief,~,~]=fieldpoints(xyb,topology,k,ChiefPoint);
    
    newA = [A;Apchief];
    newP = [inc_pressure;pIch];

    ps=[A;Apchief]\(-2*pi*[inc_pressure;pIch]);
    % Analytical solution
    pAna=cylscat(k,Rc,xyb(:,1:2),150);
    
    CylError2(nn) = mean(abs(ps-pAna.')./abs(pAna.'));
end


% POSTPROCESSING
%%
% Convergence plot
figure;
loglog(NelT*4,CylError2,'-ro'); hold on
loglog(NelT*4,CylError,'-bo');
title(['Convergence of BEM solution - k: ' num2str(k)])
xlabel('Number of elements'); ylabel('Relative error')
legend("With chief", "without chief")
grid


%% observe interior problem
