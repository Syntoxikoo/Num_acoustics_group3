% This file calculates the coefficient matrices to be used for solving. 
% It stores a set of matrices per frequency in separate files.


clear, tic

VERSION='VER3i'; % choose the calculated microphone version to process

% Input:
f=200*2.^([(0:16)/3 (17*3:30*3)/9]); % 1/3 octaves from 200 Hz to 200 kHz

m=0;
noloss=0; % if set to 1, also performs lossless calculation

% Constants
[rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu]=VTconst(f,-1);

% Load mesh and physical parameters
el_wl=ceil(6*max(f)/c); % Mesh density as a function of the frequency
eval(['[rzb38,topology38,rzline38,IndMb38,IndBP38,Mmb38,M38,N38,Tm38,mu_m38,' ...
       'rzb39,topology39,rzline39,IndMb39,IndBP39,Mmb39,M39,N39,Tm39,mu_m39,Rad,BPrad,'...
       'rzbEXT,topologyEXT,Mext,Next,IndMbEXT,MmbEXT]=geom_BK4938_9_' VERSION '(el_wl);'])

% Calculations
if noloss==1, CondNrsNL=[]; end     % The lossless pressure for each frequency

CondNrs=[]; % Condition numbers

for ii=1:length(f)
    disp(['Calculating f= ' num2str(f(ii)) ' Hz'])
    
    if noloss==1
        % Losless calculation
        [A38,B38,CConst38]=BEMEquat0(rzb38,topology38,kp(ii),m);
        CondNrsNL38=[CondNrsNL; cond(A38) cond(B38)];
        
        [A39,B39,CConst39]=BEMEquat0(rzb39,topology39,kp(ii),m);
        CondNrsNL39=[CondNrsNL; cond(A39) cond(B39)];
    end
    
     % Calculate coefficient matrices
    [Aa38,Ba38,CConst38a]=BEMEquat0(rzb38,topology38,ka(ii),m);
    [Ah38,Bh38,CConst38h]=BEMEquat0(rzb38,topology38,kh(ii),m);
    [Av38,Bv38,CConst38v]=BEMEquat0(rzb38,topology38,kv(ii),m);
%    CondNrs38=[CondNrs; cond(Aa38) cond(Ba38) cond(Ah38) cond(Bh38) cond(Av38) cond(Bv38)];

     save(['VT_BK4938-9_f' num2str(floor(f(ii)))])
%     save(['VT_BK4938-9_f' num2str(floor(f(ii)))],'-append')
% 
%     [Aext,Bext,CConstEXT]=BEMEquat0(rzbEXT,topologyEXT,k(ii),m);
% 
%     save(['VT_BK4938-9_f' num2str(floor(f(ii)))],'Aext','Bext','rzbEXT','topologyEXT','Mext','Next','IndMbEXT','MmbEXT','-append')
    
    
end

