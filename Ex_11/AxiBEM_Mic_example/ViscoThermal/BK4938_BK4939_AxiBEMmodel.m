% Simulation of the B&K 4938 and B&K 4939 microphones including the
% interior, diaphragm and exterior domain, in axisymmetrical BEM.

% This file solves the systems of equations using the pre-calculated matrices. 

clear, tic
axipath

% External sound field excitation: 1 - Axial plane wave (free field); 2 - Uniform pressure (eq. actuator)
FreeField=2; 

% Choose microphone (4938 or 4939)
MIC='4938';

% specific parameters for the diaphragms, mics [4938 4939], used in solving process.
Tm4938=3128; Tm4939=1039; % diaphragm tension, originally 
mu_m4938=8300*6.95e-6; mu_m4939=8300*2.35e-6; % surface density, density*thickness
   
% Adjustment of tension and mass of the diaphragm
if MIC=='4938'
    Tm=Tm4938;
    mu_m=mu_m4938;
else
    Tm=Tm4939;
    mu_m=mu_m4939;
end

% Excitation:
pinc=1; % EXCITATION: external (constant) incident pressure
f=200*2.^([(0:16)/3 (17*3:30*3)/9]); % 1/3 octaves from 200 Hz to 200 kHz with more frequencies
m=0; % Fully axisymmetrical

% Constants
[rho,c,kp,ka,kh,kv,tau_a,tau_h,phi_a,phi_h,eta,mu]=VTconst(f,-1);
k=2*pi*f/c;
el_wl=ceil(6*max(f)/c); % Mesh density as a function of the frequency


% INTERIOR DOMAIN (68 elements)
Hgh=19.8e-6;      % Gap thickness
Rad=2e-3;         % Diaphragm radius
BPrad=1.75e-3;    % Back plate radius
segments=[0 Hgh Rad Hgh 24 0 el_wl;...
          Rad Hgh Rad -0.75e-3 4 0 el_wl;...
          Rad -0.75e-3 1.8e-3 -0.75e-3 2 0 el_wl;...
          1.8e-3 -0.75e-3 1.8e-3 -1.5e-3 4 0 el_wl;...
          1.8e-3 -1.5e-3 1.3e-3 -1.5e-3 3 0 el_wl;...
          1.3e-3 -1.5e-3 1.3e-3 -0.5e-3 4 0 el_wl;...
          1.3e-3 -0.5e-3 BPrad -0.5e-3 3 0 el_wl;...
          BPrad -0.5e-3 BPrad 0 3 0 el_wl;...
          BPrad 0 0 0 21 0 el_wl];
[rzb,topology,rzline]=nodegen(segments,'n'); %hh38=gcf;
M=size(rzb,1);N=size(topology,1);
% indicate interior domain
rzb(:,end)=-rzb(:,end);
topology(:,end)=-topology(:,end);
% find nodes on the membrane:
IndMb=find(rzb(:,2)>Hgh-eps & rzb(:,2)<Hgh+eps & rzb(:,1)<=Rad*(1+eps));
Mmb=length(IndMb);
IndBP=find(rzb(:,2)>Hgh-eps & rzb(:,2)<Hgh+eps & rzb(:,1)<=BPrad*(1+eps));



% EXTERIOR DOMAIN:
Rext=3.175e-3;   % Exterior radius of the microphone
Lext= 50e-3;     % Length of the microphone+amplifier body
ElMemBckp=21;    % Number of elements over the membrane and backplate 
segments49EXT=[0 Hgh BPrad Hgh ElMemBckp 0 el_wl;...
          BPrad Hgh Rad Hgh round(ElMemBckp*(Rad-BPrad)/BPrad) 0 el_wl;...
          Rad Hgh Rext Hgh round(ElMemBckp*(Rext-Rad)/BPrad) 0 el_wl;...
          Rext Hgh Rext -Lext 100 0 el_wl;...
          Rext -Lext 0 -Lext-Rext 10 Rext el_wl];
[rzbEXT,topologyEXT]=nodegen(segments49EXT,'n'); %hh38=gcf;
Mext=size(rzbEXT,1);Next=size(topologyEXT,1);
% find nodes on the membrane (EXT):
IndMbEXT=find(rzbEXT(:,2)>Hgh-eps & rzbEXT(:,2)<Hgh+eps & rzbEXT(:,1)<=Rad*(1+eps));
MmbEXT=length(IndMbEXT);

[nvect38,tvect38]=normals(rzb,topology,'y');

% Constant matrices used in the calculation, a function of normal and tangential vectors:
NT1_38=nvect38(:,1)*nvect38(:,1)'+nvect38(:,2)*nvect38(:,2)'; NT2_38=nvect38(:,1)*tvect38(:,1)'+nvect38(:,2)*tvect38(:,2)';
NT3_38=tvect38(:,1)*nvect38(:,1)'+tvect38(:,2)*nvect38(:,2)'; NT4_38=tvect38(:,1)*tvect38(:,1)'+tvect38(:,2)*tvect38(:,2)';

% Gradient and laplacian matrices including the whole generator and varying spacings. Based on the secant formula applied twice.
[D0_38,Dl_38] = genderiv3(rzb,rzline);


% FEM MEMBRANE MODEL
elemem=10; % number of FEM membrane elements
nmem=2*elemem+1; % number of FEM membrane nodes
nmemBP=2*round(BPrad/Rad*(elemem))+1;
tmp=linspace(BPrad,Rad,nmem-nmemBP+1);
rmem=[linspace(0,BPrad,nmemBP) tmp(2:end)]; % membrane node mesh
topomem=[(1:2:nmem-2)' (2:2:nmem-1)' (3:2:nmem)'];  % membrane quadratic elements, for coupling interpolation
%[Am38,Bm38,rhs38]=FEMmemLIN(rmem); % Linear FEM version
[Am38,Bm38,rhs38,rnode,rtopo]=FEMmemQUAD(rmem(1:2:end)); % Quadratic FEM version. rnode and rtopo should be the same as rmem and topomem.
IndVent=find(rzb(:,1)>Rad-eps & rzb(:,2)<-0.9e-4 & rzb(:,2)>-3e-4);IndVent=IndVent(1); 


% Coupling matrices and vectors
kkk=1e-2; % Factor used on the membrane equation to reduce condition number
rhs38(end,:)=zeros(1,nmem); % Impose boundary condition, no movement at the rim (1/2)
URS38=zeros(nmem,M); % Upper Right 
for mm=1:nmem % coupling matrix to extract pressure at membrane nodes from the BEM mesh values
    % find BEM element the FEM node belongs to
    elm=find(rzb(topology(:,1),1)<=rmem(mm) & rzb(topology(:,end-1),1)>=rmem(mm) & rzb(topology(:,1),2)==rzb(1,2));elm=elm(1);
    rloc=(rmem(mm)-rzb(topology(elm,1),1))/(rzb(topology(elm,end-1),1)-rzb(topology(elm,1),1))*2-1; % local coordinate of the FEM node
    for ss=1:nmem
        URS38(ss,topology(elm,1:end-1))=URS38(ss,topology(elm,1:end-1))+rhs38(ss,mm)*[0.5*rloc.*(rloc-1) 1-rloc.^2 0.5*rloc.*(rloc+1)]/Tm/kkk;
    end
end
LLS38=zeros(M,nmem);  % Lower left
for mm=1:Mmb % coupling matrix to extract velocities at membrane nodes from the FEM mesh displacement values
    % find FEM element the BEM node belongs to
    elm=find(rmem(topomem(:,1))<=rzb(mm,1) & rmem(topomem(:,3))>=rzb(mm,1));elm=elm(1);
    rloc=(rzb(mm,1)-rmem(topomem(elm,1)))/(rmem(topomem(elm,3))-rmem(topomem(elm,1)))*2-1; % local coordinate of the BEM node
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
p38UE=zeros(nmem+M+Mext,1);
p38UE(1:nmem,1)=rhs38*ones(nmem,1)*pinc/Tm/kkk; p38UE(nmem,1)=0; % Boundary condition


% Calculations
pas38=zeros(M,length(f));     % The acoustic pressure for each frequency
phs38=zeros(M,length(f));     % The thermal pressure for each frequency
vtot38=zeros(M,length(f));    % The calculated total normal velocity on the boundary
pext=zeros(Mext,length(f));     % Exterior pressure
wmem=zeros(nmem,length(f));     % Membrane displacement
wm38=zeros(length(f),2);        % Averaged membrane displacement

% Path to pre-calculated matrices
file='.\MATfiles\VT_BK4938-9_f_VER3d\VT_BK4938-9_f';

for ii=[1:length(f)]
    % Load coefficient matrices from file
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
        p38=[zeros(nmem+M,1);Bext\pIext];
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
    
    Yv(ii)=0.1; Yvect=zeros(M,1);Yvect(IndVent,1)=Yv(ii); 
    Aa38=(Aa38+j*kp(ii)*rho*c*Ba38*diag(Yvect));
    
    % Include thermal contribution
    ABah=phi_a(ii)*(Ba38\Aa38) - phi_h(ii)*(Bh38\Ah38)*tau_a(ii)/tau_h(ii);
    
    % Include viscothermal contribution
    rho_i_0=find(abs(rzb(:,1))<1e-10); % the values where rho=0 are changed to some very small number to avoid 1/r singularity in the divergence
    rho0=rzb(:,1);
    rho0(rho_i_0)=1e-10; %%%%% This value should not be too small (instability of the system) nor too high (unrealistic geometry). Work it out as a function of the other rho values
    
    % New version of the previous corrected vectorial relationship. This one has been checked for equivalence in a dummy test 
    ABaht= inv((rho0*ones(1,M)).*NT1_38.*(Bv38\Av38) + diag(nvect38(:,1))) * C1 * ...
        (rho0*ones(1,M).*(NT2_38.*(Bv38\Av38)*D0_38) + (rho0*ones(1,M).*Dl_38) + ((tvect38(:,1)*ones(1,M)).*D0_38)); % expression with normal and tangential vectors (v2)

    CoefM=[ ULS38 URS38*(1-tau_a(ii)/tau_h(ii)) URext ;...
            -j*2*pi*f(ii)*LLS38 -(ABah+ABaht) zeros(M,Mext);...
            -j*2*pi*f(ii)*LLext zeros(Mext,M) BAext];
    BB=eye(M+nmem+Mext);
    
    % Solve for acoustic pressure
    warning off;
    pa=CoefM\(BB*p38);  % Solve the coupled system, viscous and thermal loss
    warning on;
    
    % Collect acoustic pressures for all frequencies
    pas38(:,ii)=pa(1+nmem:nmem+M);  % Viscous and thermal loss
    phs38(:,ii)=-pa(1+nmem:nmem+M)*tau_a(ii)/tau_h(ii);  % Viscous and thermal loss
    
    pext(:,ii)=pa(1+nmem+M:nmem+M+Mext);
    wmem(:,ii)=pa(1:nmem);
    vtot38(:,ii)=(ABah+ABaht)*pa(1+nmem:nmem+M);

    % Integrate displacement along radius (up to backplate) and circumference (2*pi)and divide by area over the backplate:
    RadMEM=rmem; eMEM=wmem(:,ii);
    wm38(ii,2) = trapz(RadMEM(1:nmemBP)',2*pi*RadMEM(1:nmemBP)'.*eMEM(1:nmemBP))/(pi*BPrad^2) *200/Hgh; 
    % ...and multiply by polarization voltage and divide by gap thickness (cond. mic. equations). The result is the output voltage: sensitivity.

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


figure;
if FreeField==1, Exc='axial plane wave'; else, Exc='uniform pressure';end
subplot(2,1,1)
semilogx(Mic_4938_response_f,Mic_4938_response_max+wm38dB(1),'-b',Mic_4938_response_f,Mic_4938_response_min+wm38dB(1),'-g',Mic_4938_response_f,Mic_4938_response_mean+wm38dB(1),'-k',f,wm38dB,'-rx');grid
title(['B&K ' MIC ', ' Exc ', gap: ' num2str(rzb(1,2)*1e6) ' \mum, Tm: ' num2str(Tm) ' N/m, \rho_{s}: ' num2str(mu_m*1000) ' g/m^{2}']);
legend('Actuator B&K4938 MAX','Actuator B&K4938 MIN','Actuator B&K4938 MEAN','BEM','Location','southwest')
axis([f(1) f(end) -23+wm38dB(1) 10+wm38dB(1)])
xlabel('Frequency, Hz');ylabel('Sensitivity, dB')
subplot(2,1,2)
semilogx(f,unwrap(angle(wm38(:,2)))*180/pi+180,'-rx');grid
title(['B&K ' MIC ', ' Exc ', gap: ' num2str(rzb(1,2)*1e6) ' \mum, Tm: ' num2str(Tm) ' N/m, \rho_{s}: ' num2str(mu_m*1000) ' g/m^{2}']);
legend('BEM (phase)','Location','southwest')
xxyy=axis;
axis([f(1) f(end) xxyy(3:4)])
xlabel('Frequency, Hz');ylabel('Phase, deg.')

 
% Pressure and membrane displacement
rangef=1:5:length(f);
%rangef=[1 18 33 48 57];
legendstr=cell(1,length(rangef));
for jj=1:length(rangef)
    legendstr{jj}=[num2str(floor(f(rangef(jj)))) 'Hz'];
end

figure;
subplot(2,1,1)
plot(rzline(1:Mmb,1)*(1-1e-10),abs(LLS38(1:Mmb,:)*wmem(:,rangef)),'-o');grid
%plot(rmem,abs(wmem(:,rangef)),'-o');grid
legend(legendstr)
title(['B&K ' MIC ', ' Exc ' excitation, membrane displacement'])
xlabel('Radial position, m');ylabel('Displacement, m')
subplot(2,1,2)
plot(rzline(1:Mmb,1),abs(pas38(1:Mmb,rangef)+phs38(1:Mmb,rangef)),'-x');grid  
legend(legendstr)
title(['B&K ' MIC ', ' Exc ' excitation, pressure behind the diaphragm'])
xlabel('Generator nodes');ylabel('Pressure, Pa')

