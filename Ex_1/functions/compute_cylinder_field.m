function [p_tot,p_i,p_s] = compute_cylinder_field(x,y,a,f,k,p0,M,coords,t, center,surf) % Compute the pressure field around an infinite cylinder for an incoming plane wave. 
% Input:
%   x: first coordinate (could be an int,a vector or a matrix)
%   y: second coordinate (could be an int,a vector or a matrix)
%   a: radius of the cylinder
%   f: frequency of the incoming wave
%   p0: amplitude of the incoming wave
%   M: number of terms in the series expansion
%   coords: coordinates system (cartesian or polar) of the input coordinates
%   t: time (default is 0) could be a scalar or a vector
%   center: center of the cylinder
% Output:
%   p_i: pressure field inside the cylinder
%   p_s: pressure field outside the cylinder
% -------------------------------------------------------------------------
% Example:
% ## Cartesian field
%   f=100;
%   k=0;
%   p0=1;
%   M=1000;
%   coords = 'cartesian';
%   t=0;
%   center = [0,0];
%   [p_i,p_s] = compute_cylinder_field(X,Y,a,f,k,p0,M,coords,t, center);
%   p_tot = p_s + p_i;
%   surf(X,Y,abs(p_tot));
%   ------------------------------
% ## polar field
%   [X_1,Y_1] = cart2pol(X,Y);
%   f=100;
%   k=0;
%   p0=1;
%   M=1000;
%   coords = "polar";
%   t=0;
%   cyl_ctr = [1,-2];
%   [p_i,p_s] = compute_cylinder_field(X_1,Y_1,a,f,k,p0,M,coords,t, cyl_ctr);
%   p_tot = p_s + p_i;
%   surf(X,Y,abs(p_tot));

arguments
    x 
    y 
    a = 1
    f = 0
    k = 0
    p0 = 1
    M = 1000
    coords = "cartesian"
    t = 0
    center = [0,0]
    surf = false
end
if lower(coords) == "cartesian"
    x = x-center(1);
    y = y-center(2);
    [the,r] = cart2pol(x,y);
elseif lower(coords) == "polar" 
    if center(1) == 0 && center(2) == 0 
        the = x;
        r = y;
    else % if the center is not at the origin, perform a translation
        [x,y] = pol2cart(x,y);
        x = x - center(1);
        y = y - center(2);
        [the,r] = cart2pol(x,y);
    end
end
rho = 1.21;
c0 =  343;
if k == 0
    k = 2*pi*f/c0;
elseif f == 0
    f = k*c0/(2*pi);
end

% Computing pressure field for the incoming wave
secbess = 0;
for ii = 1:M
    secbess = secbess + (1j^(ii) *cos(ii * the).* besselj(ii, k*r)); % Computing the series expansion of the Bessel function
end
p_i = p0 .* (besselj(0,k*r)+ 2 * secbess); % Computing full insident pressure field
if length(t)>1
    
    p_tmp = zeros(size(p_i,1),size(p_i,2),length(t));
    for ii = 1: length(t)
        p_tmp(:,:,ii) = p_i.* exp(-1j*2*pi*f*t(ii));
    end
    p_i = p_tmp;
end
% Computing pressure field for the scattered wave
p_s = 0;
for ii = 0:M
    if ii < 1
        eps_m =1;
        gamma_m = atan(-besselj(1,k*a)/bessely(1,k*a)); % Computing the phase angle
    else
        eps_m = 2;
        gamma_m = atan((besselj(ii-1,k*a)-besselj(ii+1,k*a))/(bessely(ii+1, k*a)- bessely(ii-1,k*a))); 
    end
    
    A_m = -eps_m*p0* 1j^(ii+1)*exp(-1j *gamma_m) * sin(gamma_m); % Computing the amplitude of the scattered wave
    if abs(A_m) <  1e-300
        break
    end
    if isnan(A_m)
        disp("error for m =" + ii)
        break
    end

    p_s = p_s+ (A_m *cos(ii * the).*( besselj(ii, k*r)+1j * bessely(ii,k*r)));
end
if length(t)>1
    p_tmp = zeros(size(p_s,1),size(p_s,2),length(t));
    for ii = 1: length(t)
        p_tmp(:,:,ii) = p_s.* exp(-1j*2*pi*f*t(ii));
    end
    p_s = p_tmp;
end
if surf == false
    mask = find(r <= a);
    p_i(mask) = NaN;
    p_s(mask) = NaN;
end
p_tot = p_s + p_i;
end