function IPsing=singgausstri(n,nodn,varargin)
% IPsing=singgausstri(n,nodn,see)
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

% Round n to the nearest even integer
n=ceil(n/2)*2;

if any(nodn<=3)
   n_ang=n; n_rad=n;
   [bp_ang,wf_ang]=gaussrule(n_ang);
   [bp_rad,wf_rad]=gaussrule(n_rad);
   
   a=0;b=pi/2;
   theta=(a+b)/2+(b-a)/2*bp_ang;
   ang_weight=(b-a)/2*wf_ang;

   aa=0;bb=1./(cos(theta)+sin(theta));
   rad=ones(n_rad,1)*(aa+bb)'/2+bp_rad*(bb-aa)'/2; % positions along the radii
   rad_weight=(wf_rad*(ang_weight.*(bb-aa)/2)').*rad; % weights along the radii
   % (rad is the jacobean of the polar transformation)
   
   % Transform back from polar to local (x,y) coordinate system
   IPa=[rad.*(ones(n_rad,1)*cos(theta)') rad.*(ones(n_rad,1)*sin(theta)') rad_weight];
        
   IPa=[IPa(1:end/3)' IPa(end/3+1:2*end/3)' IPa(2*end/3+1:end)'];
end
   
if any(nodn>=4)
   n_ang=n/2; n_rad=n;
   [bp_ang,wf_ang]=gaussrule(n_ang);
   [bp_rad,wf_rad]=gaussrule(n_rad);
   
   % First do the triangle from theta = 0 to atan(2)
   a=0;b=atan(2);
   theta=(a+b)/2+(b-a)/2*bp_ang;
   ang_weight=(b-a)/2*wf_ang;
   aa=0;bb=0.5./cos(theta);
   rad=ones(n_rad,1)*(aa+bb)'/2+bp_rad*(bb-aa)'/2; % positions along the radii
   rad_weight=(wf_rad*(ang_weight.*(bb-aa)/2)').*rad; % weights along the radii
   % (rad is the jacobean of the polar transformation)
   
   % Transform back from polar to local (x,y) coordinate system
   xx=rad.*(ones(n_rad,1)*cos(pi-theta)')+0.5;
   yy=rad.*(ones(n_rad,1)*sin(pi-theta)');
   IPb=[xx(1:end)' yy(1:end)' rad_weight(1:end)'];
   
   % Then do the triangle from atan(2) to pi
   a=atan(2);b=pi;
   theta=(a+b)/2+(b-a)/2*bp_ang;
   ang_weight=(b-a)/2*wf_ang;
   aa=0;bb=sin(pi/4)./sin(theta-pi/4)*0.5;
   rad=ones(n_rad,1)*(aa+bb)'/2+bp_rad*(bb-aa)'/2; % positions along the radii
   rad_weight=(wf_rad*(ang_weight.*(bb-aa)/2)').*rad; % weights along the radii
   % (rad is the jacobean of the polar transformation)
   
   % Transform back from polar to local (x,y) coordinate system
   xx=rad.*(ones(n_rad,1)*cos(pi-theta)')+0.5;
   yy=rad.*(ones(n_rad,1)*sin(pi-theta)');
   IPb=[IPb ;xx(1:end)' yy(1:end)' rad_weight(1:end)'];
end

for nn=1:length(nodn)
    switch nodn(nn)
        case 1
            IPsing(:,:,nn)=[1-IPa(:,1)-IPa(:,2) IPa(:,1) IPa(:,3)];
        case 2
            IPsing(:,:,nn)=[IPa(:,2) 1-IPa(:,1)-IPa(:,2) IPa(:,3)];
        case 3
            IPsing(:,:,nn)=IPa;
        case 4
            IPsing(:,:,nn)=[1-IPb(:,1)-IPb(:,2) IPb(:,1) IPb(:,3)];
        case 5
            IPsing(:,:,nn)=[IPb(:,2) 1-IPb(:,1)-IPb(:,2) IPb(:,3)];
        case 6
            IPsing(:,:,nn)=IPb;
    end
end

% To see the result:
if nargin>2 & varargin{1}=='y'
    figure;
    for nn=1:length(nodn)
        subplot(floor(((length(nodn))-1)/3)+1,3,nn)
        plot3(IPsing(:,1,nn),IPsing(:,2,nn),IPsing(:,3,nn),'.');
        %    axis([0 1 0 1])
        view(0,90);
    end
end

end
