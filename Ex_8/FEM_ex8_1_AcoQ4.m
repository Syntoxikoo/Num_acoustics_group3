function [Ke,Me]=FEM_ex8_1_AcoQ4(ex,ey,ep)

%% -------------------------------------------------------------
% PURPOSE
%  Compute element stiffness and mass  
%  matrix for 4 node isoparametric acoustic element
%
% INPUT:  ex = [x1 x2 x3 x4]   element coordinates
%         ey = [y1 y2 y3 y4]                           
%         ep = [c rho0 ir]   speed of sound, density, integration rule
%
% OUTPUT: Ke :  element 'stiffness' matrix (4 x 4)
%         Me :  element mass matrix (4 x 4)
% Refer from CalFEM by G. Sandberg, P.-A. Wernberg and P. Davidsson
%% -------------------------------------------------------------

c = ep(1); 
rho0 = ep(2); 
ir=ep(3); 
ngp = ir*ir;

%% Gauss points and its weights 
if ir==1
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
elseif ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
elseif ir==3
    g1=0.774596699241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;
    gp(:,1)=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2)=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];
else
    disp('Used number of integration points not implemented');
    return
end
wp=w(:,1).*w(:,2);

%% Shape function and deriviatve of shape function with respect to u,v
  u=gp(:,1);  v=gp(:,2);  r2=ngp*2;

  N(:,1)=(1-u).*(1-v)/4;  
  N(:,2)=(1+u).*(1-v)/4;
  N(:,3)=(1+u).*(1+v)/4; 
  N(:,4)=(1-u).*(1+v)/4;

  dNr(1:2:r2,1)=-(1-v)/4; 
  dNr(1:2:r2,2)= (1-v)/4;
  dNr(1:2:r2,3)= (1+v)/4;    
  dNr(1:2:r2,4)=-(1+v)/4;
  
  dNr(2:2:r2+1,1)=-(1-u)/4;  
  dNr(2:2:r2+1,2)=-(1+u)/4;
  dNr(2:2:r2+1,3)= (1+u)/4;  
  dNr(2:2:r2+1,4)= (1-u)/4;


  Ke1=zeros(4,4); 
  Me1=zeros(4,4);  
  
  %% Jacobian Matrix; 
  JT=dNr*[ex;ey]';


  for i=1:ngp
    idx = [ 2*i-1; 2*i ]; 
    detJ = JT(idx(1),1)*JT(idx(2),2) - JT(idx(2),1)*JT(idx(1),2); % Determinant Jacobian
    B    = inv(JT(idx,:)) * dNr(idx,:)  ; % B matrix
    % Notes: Weights are safed in the vector wp
    Me1  = Me1 + N(i,:).'*N(i,:)*detJ*wp(i); % Element mass matrix
    Ke1  = Ke1 + B.'*B *detJ*wp(i); % Element Stiffness matrix
  end
  
  Me=Me1;
  Ke=Ke1*c^2; 

end

