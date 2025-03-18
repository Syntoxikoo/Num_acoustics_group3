a =1; %radius
the = linspace(0,-2*pi, 360);
% r=1;
k= 1:5;
p0=1;
f=0;
M=1000;
coords = 'polar';
t=0;
cyl_ctr = [0,0];
p_i = zeros(length(the),length(k));
p_s = zeros(length(the),length(k));
for ii = 1:length(k)
    r =a*10;
    [p_i(:,ii),p_s(:,ii)] = compute_cylinder_field(the,r,a,f,k(ii),p0,M,coords,t, cyl_ctr);
end
