% BK4938

% Includes EXTERIOR domain

% This file solves the systems of equations using the pre-calculated matrices. 
% It obtains the pressures and velocities to be processed later to show the results.

% MinMEM: This version reduces the size of the membrane FEM matrices, using
% fewer membrane points. The condition number is also reduced.
% Includes two implementations of a vent, as a pressure release node and as
% a prescribed impedance node. The impedance is defined in a special
% function using formulas in the literature.
% The EXT version couples interior, membrane and external body. The
% exterior and interior BEM do not need matching meshes, so an older
% version is used for the exterior.

clear, tic

% Geometry version and coefficient matrices to be loaded:
s_ver='';  % sub-version
 VERSION='VER3d'; % Hgh=19.8e-6, 68 elements, sharp corners    %  
% VERSION='VER3c'; % Hgh=19.8e-6, 108 el, rounded corners    %
% VERSION='VER3e'; % Hgh=19.65e-6, 68 elements, sharp corners  
% VERSION='VER3f'; % Hgh=18e-6, 68 elements, sharp corners
% VERSION='VER3g'; % Hgh=19.5e-6, 68 elements, sharp corners     
% VERSION='VER3h'; % Hgh=19.6e-6, 68 elements, sharp corners
% VERSION='VER3i'; % Hgh=19.55e-6, 68 elements, sharp corners    
% VERSION='VER4b'; % Hgh=   20e-6, 153 el, rounded corners (more gap elements)
% VERSION='VER4c'; % Hgh=19.8e-6, 153 el, rounded corners (more gap elements)  %

% Input:
pinc=1; % EXCITATION: external (constant) incident pressure
f=200*2.^([(0:16)/3 (17*3:30*3)/9]); % 1/3 octaves from 200 Hz to 200 kHz with more frequencies

FreeField=0; % If 1, a free field excitation (axial plane wave) should be used. In this file the external excitation is prescribed as in an actuator.
m=0;

% Constants
[rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu]=VTconst(f,-1);
k=2*pi*f/c;

% Load mesh and physical parameters
el_wl=ceil(6*max(f)/c); % Mesh density as a function of the frequency
file=['.\MATfiles\VT_BK4938-9_f_' VERSION s_ver '\VT_BK4938-9_f'];
eval(['[rzb38,topology38,rzline38,IndMb38,IndBP38,Mmb38,M38,N38,Tm38,mu_m38,' ...
       'rzb39,topology39,rzline39,IndMb39,IndBP39,Mmb39,M39,N39,Tm39,mu_m39,Rad,BPrad,' ...
          'rzbEXT,topologyEXT,Mext,Next,IndMbEXT,MmbEXT,Hgh]=geom_BK4938_9_' VERSION '(el_wl);'])

% Load parameters of the exterior geometry
topologyEXT=[]; % For some unknown reason, Matlab wants this variable to be pre-defined
eval(['[dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,dum,'...
       'rzbEXT,topologyEXT,Mext,Next,IndMbEXT,MmbEXT]=geom_BK4938_9_VER3(el_wl,''a'');'])


% specific parameters for the diaphragms, mics [4938 4939], used in solving process.
Tm38=3128; Tm39=1039;
mu_m38=8300*6.95e-6; mu_m39=8300*2.35e-6;% surface density, density*thickness 8300*[6.95e-6 2.35e-6]   
   
% Adjustment of tension and mass of the diaphragm
 Tm=Tm38*1.1;
 mu_m=mu_m38;

[nvect38,tvect38]=normals(rzb38,topology38,'y');

% Constant matrices used in the calculation, a function of normal and tangential vectors:
NT1_38=nvect38(:,1)*nvect38(:,1)'+nvect38(:,2)*nvect38(:,2)'; NT2_38=nvect38(:,1)*tvect38(:,1)'+nvect38(:,2)*tvect38(:,2)';
NT3_38=tvect38(:,1)*nvect38(:,1)'+tvect38(:,2)*nvect38(:,2)'; NT4_38=tvect38(:,1)*tvect38(:,1)'+tvect38(:,2)*tvect38(:,2)';

% Gradient and laplacian matrices including the whole generator and varying spacings. Based on the secant formula applied twice.
[D0_38,Dl_38] = genderiv3(rzb38,rzline38);


% FEM MEMBRANE MODEL
elemem=10; % number of FEM membrane elements
nmem=2*elemem+1; % number of FEM membrane nodes
nmemBP=2*round(BPrad/Rad*(elemem))+1;
tmp=linspace(BPrad,Rad,nmem-nmemBP+1);
rmem=[linspace(0,BPrad,nmemBP) tmp(2:end)]; % membrane node mesh
topomem=[(1:2:nmem-2)' (2:2:nmem-1)' (3:2:nmem)'];  % membrane quadratic elements, for coupling interpolation
%[Am38,Bm38,rhs38]=FEMmemLIN(rmem); % Linear FEM version
[Am38,Bm38,rhs38,rnode,rtopo]=FEMmemQUAD(rmem(1:2:end)); % Quadratic FEM version. rnode and rtopo should be the same as rmem and topomem.
IndVent=find(rzb38(:,1)>Rad-eps & rzb38(:,2)<-0.9e-4 & rzb38(:,2)>-3e-4);IndVent=IndVent(1); 


% Coupling matrices and vectors (4938);
kkk=1e-2; % Factor used on the membrane equation to reduce condition number
rhs38(end,:)=zeros(1,nmem); % Impose boundary condition, no movement at the rim (1/2)
URS38=zeros(nmem,M38); % Upper Right 
for mm=1:nmem % coupling matrix to extract pressure at membrane nodes from the BEM mesh values
    % find BEM element the FEM node belongs to
    elm=find(rzb38(topology38(:,1),1)<=rmem(mm) & rzb38(topology38(:,end-1),1)>=rmem(mm) & rzb38(topology38(:,1),2)==rzb38(1,2));elm=elm(1);
    rloc=(rmem(mm)-rzb38(topology38(elm,1),1))/(rzb38(topology38(elm,end-1),1)-rzb38(topology38(elm,1),1))*2-1; % local coordinate of the FEM node
    for ss=1:nmem
        URS38(ss,topology38(elm,1:end-1))=URS38(ss,topology38(elm,1:end-1))+rhs38(ss,mm)*[0.5*rloc.*(rloc-1) 1-rloc.^2 0.5*rloc.*(rloc+1)]/Tm/kkk;
    end
end
LLS38=zeros(M38,nmem);  % Lower left
for mm=1:Mmb38 % coupling matrix to extract velocities at membrane nodes from the FEM mesh displacement values
    % find FEM element the BEM node belongs to
    elm=find(rmem(topomem(:,1))<=rzb38(mm,1) & rmem(topomem(:,3))>=rzb38(mm,1));elm=elm(1);
    rloc=(rzb38(mm,1)-rmem(topomem(elm,1)))/(rmem(topomem(elm,3))-rmem(topomem(elm,1)))*2-1; % local coordinate of the BEM node
    LLS38(mm,topomem(elm,1:3))=[0.5*rloc.*(rloc-1) 1-rloc.^2 0.5*rloc.*(rloc+1)];
end


% Coupling matrices and vectors (Exterior);
URext=zeros(nmem,Mext); % Upper Right 
for mm=1:nmem % coupling matrix to extract pressure at membrane nodes from the exterior BEM mesh values
    % find BEM element the FEM node belongs to
    elm=find(rzbEXT(topologyEXT(:,1),1)<=rmem(mm) & rzbEXT(topologyEXT(:,end-1),1)>=rmem(mm) & rzbEXT(topologyEXT(:,1),2)==rzbEXT(1,2));elm=elm(1);
    rloc=(rmem(mm)-rzbEXT(topologyEXT(elm,1),1))/(rzbEXT(topologyEXT(elm,end-1),1)-rzbEXT(topologyEXT(elm,1),1))*2-1; % local coordinate of the FEM node
    for ss=1:nmem
        URext(ss,topologyEXT(elm,1:end-1))=URext(ss,topologyEXT(elm,1:end-1))+rhs38(ss,mm)*[0.5*rloc.*(rloc-1) 1-rloc.^2 0.5*rloc.*(rloc+1)]/Tm/kkk;
    end
end
LLext=zeros(Mext,nmem);  % Lower left
for mm=1:MmbEXT % coupling matrix to extract exterior velocities at membrane nodes from the FEM mesh displacement values
    % find FEM element the BEM node belongs to
    elm=find(rmem(topomem(:,1))<=rzbEXT(mm,1) & rmem(topomem(:,3))>=rzbEXT(mm,1));elm=elm(1);
    rloc=(rzbEXT(mm,1)-rmem(topomem(elm,1)))/(rmem(topomem(elm,3))-rmem(topomem(elm,1)))*2-1; % local coordinate of the BEM node
    LLext(mm,topomem(elm,1:3))=[0.5*rloc.*(rloc-1) 1-rloc.^2 0.5*rloc.*(rloc+1)];
end

% uniform excitation
p38UE=zeros(nmem+M38+Mext,1);
p38UE(1:nmem,1)=rhs38*ones(nmem,1)*pinc/Tm/kkk; p38UE(nmem,1)=0; % Boundary condition

% [Yavg38]=AVGmem([rmem(1:nmemBP)' zeros(nmemBP,1) ones(nmemBP,1)]);

% Calculations
pas38=zeros(M38,length(f));     % The acoustic pressure for each frequency
phs38=zeros(M38,length(f));     % The thermal pressure for each frequency
vtot38=zeros(M38,length(f));    % The calculated total normal velocity on the boundary
pext=zeros(Mext,length(f));     % Exterior pressure
wmem=zeros(nmem,length(f));     % Membrane displacement
wm38=zeros(length(f),2);        % Averaged membrane displacement

for ii=[1:length(f)]
    % Load coefficient matrices from file generated by 'VT_BK4938_9_CalculateMatrices.m'
    SS=load([file num2str(floor(f(ii)))],'Aa38','Ba38','Ah38','Bh38','Av38','Bv38'); % Load all versions of A and B matrices for f(ii):
    Aa38=SS.Aa38;Ba38=SS.Ba38;Ah38=SS.Ah38;Bh38=SS.Bh38;Av38=SS.Av38;Bv38=SS.Bv38;

    disp(['Processing f= ' num2str(f(ii)) ' Hz '])
    
    % Exterior domain matrices
    GG=load(['.\MATfiles\VT_BK4938-9_f_EXT\VT_BK4938-9_f' num2str(floor(f(ii)))],'Aext','Bext'); % Load Aext and Bext of the exterior for f(ii)
    Aext=GG.Aext;Bext=GG.Bext;
    Bext=-i*kp(ii)*rho*c*Bext; % B times i*k*rho*c to compute pressure instead of velocity potential
    BAext=Bext\Aext;
    
    if FreeField==1
        % Incoming wave, free field case
        pIext=incoming([0 0 pinc 0],rzbEXT,k(ii),m);
        p38=[zeros(nmem+M38,1);Bext\pIext];
    else
        p38=p38UE;
    end
    
    
    % FEM membrane 
    % membrane wavenumber, w/c: see Morse (Vibration and Sound, Ch:5, Sec 17), Robey1954, Zuckerwar 1978
    K2_38=(2*pi*f(ii)).^2*mu_m/Tm; 
    ULS38=(Am38+K2_38*Bm38)/kkk; % left hand side
    ULS38(end,:)=[zeros(1,nmem-1) 1];  % Impose boundary condition, no movement at the rim (2/2)
    ULS38(:,end)=[zeros(nmem-1,1) ; 1];  % Impose boundary condition, no movement at the rim (2/2 extra)
    
    % Matrix multiplications
    C1=phi_a(ii)-phi_h(ii)*tau_a(ii)/tau_h(ii);
    AvBv=Av38\Bv38; % At^-1*Bt
    
    Yv(ii)=0.1; Yvect=zeros(M38,1);Yvect(IndVent,1)=Yv(ii); 
    Aa38=(Aa38+j*kp(ii)*rho*c*Ba38*diag(Yvect));
    
    % Include thermal contribution
    ABah=phi_a(ii)*(Ba38\Aa38) - phi_h(ii)*(Bh38\Ah38)*tau_a(ii)/tau_h(ii);
    
    % Include viscothermal contribution
    rho_i_0=find(abs(rzb38(:,1))<1e-10); % the values where rho=0 are changed to some very small number to avoid 1/r singularity in the divergence
    rho0=rzb38(:,1);
    rho0(rho_i_0)=1e-10; %%%%% This value should not be too small (instability of the system) nor too high (unrealistic geometry). Work it out as a function of the other rho values
    
    % New version of the previous corrected vectorial relationship. This one has been checked for equivalence in a dummy test 
    ABaht= inv((rho0*ones(1,M38)).*NT1_38.*(Bv38\Av38) + diag(nvect38(:,1))) * C1 * ...
        (rho0*ones(1,M38).*(NT2_38.*(Bv38\Av38)*D0_38) + (rho0*ones(1,M38).*Dl_38) + ((tvect38(:,1)*ones(1,M38)).*D0_38)); % expression with normal and tangential vectors (v2)

    CoefM=[ ULS38 URS38*(1-tau_a(ii)/tau_h(ii)) URext ;...
            -j*2*pi*f(ii)*LLS38 -(ABah+ABaht) zeros(M38,Mext);...
            -j*2*pi*f(ii)*LLext zeros(Mext,M38) BAext];
    BB=eye(M38+nmem+Mext);
    
    % Solve for acoustic pressure
    warning off;
    pa=CoefM\(BB*p38);  % Solve the coupled system, viscous and thermal loss
    warning on;
    
    % Collect acoustic pressures for all frequencies
    pas38(:,ii)=pa(1+nmem:nmem+M38);  % Viscous and thermal loss
    phs38(:,ii)=-pa(1+nmem:nmem+M38)*tau_a(ii)/tau_h(ii);  % Viscous and thermal loss
    
    pext(:,ii)=pa(1+nmem+M38:nmem+M38+Mext);
    wmem(:,ii)=pa(1:nmem);
    vtot38(:,ii)=(ABah+ABaht)*pa(1+nmem:nmem+M38);

    % Integrate displacement along radius (up to backplate) and circumference (2*pi)and divide by area over the backplate:
    RadMEM=rmem; eMEM=wmem(:,ii);
    wm38(ii,2) = trapz(RadMEM(1:nmemBP)',2*pi*RadMEM(1:nmemBP)'.*eMEM(1:nmemBP))/(pi*BPrad^2) *200/Hgh; 
    % ...and multiply by polarization voltage and divide by gap thickness (cond. mic. equations). The result is the output voltage: sensitivity.
%     wm38(ii,2) = trapz(rmem(1:nmemBP)',2*pi*rmem(1:nmemBP)'.*wmem(1:nmemBP,ii))/(pi*BPrad^2); % Version not suitable for "parfor"
    
%     % Old version, wrong result
%     wm38(ii,2) = [Yavg38] * (pa(1:nmemBP));

end

% VERIFY RESULTS
% wm38dB=20*log10(abs(wm38(:,2))/abs(wm38(1,2)));
wm38dB=20*log10(abs(wm38(:,2)));

% Load B&K actuator measurement data
load Mic_4938_response_spread
Mic_4938_response_f=Mic_4938_response_spread(1,2:end);              % Frequency axis
Mic_4938_response_max=max(Mic_4938_response_spread(2:end,2:end));   % Maximum values
Mic_4938_response_min=min(Mic_4938_response_spread(2:end,2:end));   % Minimum values
Mic_4938_response_mean=mean(Mic_4938_response_spread(2:end,2:end)); % Average values


figure;subplot(2,1,1)
semilogx(Mic_4938_response_f,Mic_4938_response_max+wm38dB(1),'-b',Mic_4938_response_f,Mic_4938_response_min+wm38dB(1),'-g',Mic_4938_response_f,Mic_4938_response_mean+wm38dB(1),'-k',f,wm38dB,'-rx');grid
title(['MIC 4938, gap: ' num2str(rzb38(1,2)*1e6) ' \mum, Tm: ' num2str(Tm) ' N/m, \rho_{s}: ' num2str(mu_m*1000) ' g/m^{2}']);
legend('Measurement MAX','Measurement MIN','Measurement MEAN','BEM')
axis([f(1) f(end) -23+wm38dB(1) 2+wm38dB(1)])
xlabel('Frequency, Hz');ylabel('Sensitivity, dB')
subplot(2,1,2);
semilogx(Mic_4938_response_f,Mic_4938_response_max+wm38dB(1),'-b',Mic_4938_response_f,Mic_4938_response_min+wm38dB(1),'-g',Mic_4938_response_f,Mic_4938_response_mean+wm38dB(1),'-k',f,wm38dB,'-rx');grid
title('MIC 4938 (Zoom in)');legend('Measurement MAX','Measurement MIN','Measurement MEAN','BEM')
axis([4000 100000 -5+wm38dB(1) 5+wm38dB(1)])
xlabel('Frequency, Hz');ylabel('Sensitivity, dB')

 
Mmb=Mmb38;
% Pressure and membrane displacement
rangef=1:5:length(f);
%rangef=[1 18 33 48 57];
legendstr=cell(1,length(rangef));
for jj=1:length(rangef)
    legendstr{jj}=[num2str(floor(f(rangef(jj)))) 'Hz'];
end

figure;
subplot(2,1,1)
plot(rzline38(1:Mmb,1)*(1-1e-10),abs(LLS38(1:Mmb,:)*wmem(:,rangef)),'-o');grid
%plot(rmem,abs(wmem(:,rangef)),'-o');grid
legend(legendstr)
title('Membrane displacement 4938')
xlabel('Radial position, m');ylabel('Displacement, m')
subplot(2,1,2)
plot(rzline38(:,1),abs(pas38(:,rangef)+phs38(:,rangef)),'-x');grid  
legend(legendstr)
title('Pressure along the generator 4938')
xlabel('Generator nodes');ylabel('Pressure, Pa')

