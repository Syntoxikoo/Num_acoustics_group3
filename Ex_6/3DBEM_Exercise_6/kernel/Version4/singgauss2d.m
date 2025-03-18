function IPsing=singgauss2d(n,nodn)
% IPsing=singgauss2d(n)
% Produces (x,y) coordinates and weights for 2-D Gauss-Legendre quadrature
% by transforming into polar system (r,theta) where r=0 at the nodes given 
% in the vector nodn.
% The Jacobean r is worked into the weights and this formulae is well
% suited for integrating 1/R singularities at local nodes.
%
% Input parameters:
%   -n:      order of formulae in each r,theta direction
%            (if n uneven theta direction is n+1).
%   -nodn:   node number(s) with singularities. One different set of points
%            and weights is calculated for each of them.
%
% Output parameters:
%   -IPsing: n^2 x 3 x length(nodn) matrix.
%            1st column x-values; 2nd column y-values; 3rd column weigths.
%            One set in each 3rd-dimension index for each node in 'nodn'.

if any(nodn<=4)
   n_ang=ceil(n/2); n_rad=n;
   [bp_ang,wf_ang]=gaussrule(n_ang);
   [bp_rad,wf_rad]=gaussrule(n_rad);
   
   a=0;b=pi/4;
   theta=(a+b)/2+(b-a)/2*bp_ang;
   ang_weight=(b-a)/2*wf_ang;

   aa=0;bb=2./cos(theta);
   rad=ones(n_rad,1)*(aa+bb)'/2+bp_rad*(bb-aa)'/2; % positions along the radii
   rad_weight=(wf_rad*(ang_weight.*(bb-aa)/2)').*rad; % weights along the radii
   % (rad is the jacobean of the polar transformation)
   
   % Transform back from polar to local (x,y) coordinate system
   IPa=[rad.*(ones(n_rad,1)*cos(theta)')-1 rad.*(ones(n_rad,1)*sin(theta)')-1 rad_weight ;...
        % Take care of the triangle from pi/4 to pi/2
        rad.*(ones(n_rad,1)*sin(theta)')-1 rad.*(ones(n_rad,1)*cos(theta)')-1 rad_weight];
        
   IPa=[IPa(1:end/3)' IPa(end/3+1:2*end/3)' IPa(2*end/3+1:end)'];
end
   
if any(nodn>=5)
   n_ang=ceil(n/4); n_rad=n;
   [bp_ang,wf_ang]=gaussrule(n_ang);
   [bp_rad,wf_rad]=gaussrule(n_rad);
   
   % First do the triangle from theta = 0 to atan(2)
   a=0;b=atan(2);
   theta=(a+b)/2+(b-a)/2*bp_ang;
   ang_weight=(b-a)/2*wf_ang;
   aa=0;bb=1./cos(theta);
   rad=ones(n_rad,1)*(aa+bb)'/2+bp_rad*(bb-aa)'/2; % positions along the radii
   rad_weight=(wf_rad*(ang_weight.*(bb-aa)/2)').*rad; % weights along the radii
   % (rad is the jacobean of the polar transformation)
   
   % Transform back from polar to local (x,y) coordinate system
   IPb=[rad.*(ones(n_rad,1)*cos(theta)') rad.*(ones(n_rad,1)*sin(theta)')-1 rad_weight];
   
   % Then do the triangle from atan(2) to pi/2
   a=atan(2);b=pi/2;
   theta=(a+b)/2+(b-a)/2*bp_ang;
   ang_weight=(b-a)/2*wf_ang;
   aa=0;bb=2./sin(theta);
   rad=ones(n_rad,1)*(aa+bb)'/2+bp_rad*(bb-aa)'/2; % positions along the radii
   rad_weight=(wf_rad*(ang_weight.*(bb-aa)/2)').*rad; % weights along the radii
   % (rad is the jacobean of the polar transformation)
   
   % Transform back from polar to local (x,y) coordinate system
   IPb=[IPb ; rad.*(ones(n_rad,1)*cos(theta)') rad.*(ones(n_rad,1)*sin(theta)')-1 rad_weight];
   
   IPb=[IPb(1:end/3)' IPb(end/3+1:2*end/3)' IPb(2*end/3+1:end)']; % re-order
   IPb=[IPb; -IPb(:,1) IPb(:,2) IPb(:,3)]; % create the other half, from pi/2 to pi
end
   
for nn=1:length(nodn)
   if nodn(nn)<=4
      IPsing(:,:,nn)=[IPa(:,1).*sign((nodn(nn)-1.5).*(nodn(nn)-3.5)) ...
                      IPa(:,2).*sign((nodn(nn)-2.5).*(nodn(nn)-4.5)) ...
                      IPa(:,3)];
   else
      IPsing(:,:,nn)=[IPb(:,1) ...
                      IPb(:,2).*sign((nodn(nn)-5.5).*(nodn(nn)-7.5)) ...
                      IPb(:,3)];
      if nodn(nn)==6|nodn(nn)==8
         IPsing(:,:,nn)=[IPsing(:,2,nn) IPsing(:,1,nn) IPsing(:,3,nn)];
      end
   end
end

% To see the result:
%for nn=1:length(nodn)
%   figure;
%   plot3(IPsing(:,1,nn),IPsing(:,2,nn),IPsing(:,3,nn),'.');
%   view(0,90);
%end
