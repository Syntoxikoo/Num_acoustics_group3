% This file calculates the coefficient matrices to be used for solving. 
% It stores a set of matrices per frequency in separate files.
% Version made specifically for the Version 7 case


clear, tic

VERSION='VER7'; % choose the calculated microphone version to process
SubVer='n';
% Version 7: Idealized microphone with no back cavity in VCH thesis and Kampinga thesis.

% Input:
f=200*2.^([(0:16)/3 (17*3:30*3)/9]); % 1/3 octaves from 200 Hz to 200 kHz

m=0;
noloss=0; % if set to 1, also performs lossless calculation


[rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu]=VTconst(f,-1);

% Load mesh and physical parameters
el_wl=ceil(6*max(f)/c); % Mesh density as a function of the frequency
[rzb,topology,rzline,IndMb,IndSide,Mmb,M,N,Tm,mu_m,Rad,BPrad]=geom_BK4938_9_VER7(el_wl,SubVer);

% Calculations
if noloss==1, CondNrsNL=[]; end     % The lossless pressure for each frequency

CondNrs=[]; % Condition numbers

ii=1;
for ii=1:length(f)
    disp(['Calculating f= ' num2str(f(ii)) ' Hz'])
    
    if noloss==1
        % Losless calculation
        [A,B,CConst]=BEMEquat0(rzb,topology,kp(ii),m);
        CondNrsNL=[CondNrsNL; cond(A) cond(B)];
    end
    
     % Calculate coefficient matrices
    [Aa,Ba,CConsta]=BEMEquat0(rzb,topology,ka(ii),m);
    [Ah,Bh,CConsth]=BEMEquat0(rzb,topology,kh(ii),m);
    [Av,Bv,CConstv]=BEMEquat0(rzb,topology,kv(ii),m);
    CondNrs=[CondNrs; cond(Aa) cond(Ba) cond(Ah) cond(Bh) cond(Av) cond(Bv)];
    
   
    save(['VT_IdealMIC_f' num2str(floor(f(ii)))])

    
end


