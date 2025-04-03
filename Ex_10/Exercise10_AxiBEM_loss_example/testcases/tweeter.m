function [tweetersegments,membr_min,membr_mid,membr_max,tweeter_end] = tweeter(lambda,varargin);

% [tweetersegments,membr_min,membr_mid,membr_max,tweeter_end] = tweeter(lambda,Zrim);
%
% Geometrical definition of the Vifa tweeter, ready for meshing by "nodegen".
%
% Input:
%     -lambda:   wavelength.
%     -Zrim:     z value at the rim. If not specified, it is 5 mm.
%
% Output:
%     -tweetersegments:  ready for nodegen.
%     -membr_min:        [rho z] start of the membrane
%     -membr_mid:        [rho z] middle of the membrane
%     -membr_max:        [rho z] end of the membrane
%     -tweeter_end:      [rho z] end of the tweeter

% Johan Heilmann 4-2007

if nargin>1
   Zrim=varargin{1}-5e-3;
else
   Zrim=0;
end

% MIDDLESECTION OF THE TWEETER
rho1=0;                     z1=12e-3 + Zrim;   R1=8e-3;
rho2=sqrt(15)*5e-3-15e-3;   z2=6.5e-3 + Zrim;  R2=20e-3;

% MEMBRANE
rho3=5e-3;      z3=1.5e-3 + Zrim;  R3=4.8e-3;
rho4=13.05e-3;  z4= Zrim;       R4=3.34e-3;
rho5=19.5e-3;   z5= Zrim;       R5=4e-3;

% EDGE OF THE TWEETER
rho6=20.5e-3;   z6=2.1875e-3 + Zrim;   R6=-9e-3;
rho7=23.5e-3;   z7=4.0625e-3 + Zrim;   R7=46.63e-3;
rho8=34e-3;     z8=5e-3 + Zrim;        

% CONSTRUCTION OF THE TWEETER
%
% segments=[Rho_1 Z_1 Rho_2 Z_2 No._of_elements Curvature_radius Segment_admittance;
%    Rho_2 Z_2 Rho_3 Z_3 No._of_elements Curvature_radius Segment_admittance;... ];

tweetersegments=[rho1 z1 rho2 z2 5 R1 6/lambda;
    rho2 z2 rho3 z3 4 R2 6/lambda;
    rho3 z3 rho4 z4 6 R3 6/lambda;
    rho4 z4 rho5 z5 6 R4 6/lambda;
    rho5 z5 rho6 z6 3 R5 6/lambda;
    rho6 z6 rho7 z7 3 R6 6/lambda;
    rho7 z7 rho8 z8 3 R7 6/lambda];

 
% SPECIAL POINTS OF INTEREST 
tweeter_end=[rho8 z8];  % End of the tweeter
membr_min=[rho3 z3];    % Start of the membrane
membr_mid=[rho4 z4];    % Middle of the membrane
membr_max=[rho5 z5];    % End of the membrane



% membrane1segments=[rho3 z3 rho4 z4 max(6,ceil(pi*R3*6/lambda)) R3 0];
% membrane2segments=[rho4 z4 rho5 z5 max(6,ceil(pi*R4*6/lambda)) R4 0];
% [rzb1,topology]=nodegen(membrane1segments,'y');
% [rzb2,topology]=nodegen(membrane2segments,'y');
% 
% M1=size(rzb1,1);
% N1=size(topology1,1);
% M2=size(rzb2,1);
% N2=size(topology2,1);
% 
% nn1=find(rzb1(:,1)>=membr_min(1) & rzb1(:,1)<=membr_mid(1));
% nn2=find(rzb2(:,1)>=membr_mid(1) & rzb2(:,1)<=membr_max(1));
% 
% 
% zz=inline('a*rho+b');
% left=zz(a1,b1,);
% right=zz(a2,b2,);
% 
% % vs=zeros(M,1); vs(nn)=ones(length(nn),1)*u0; % !!!!!!!!!!!!!!!
% 
% 
% 
