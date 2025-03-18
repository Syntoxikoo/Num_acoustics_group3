x = linspace(-10,10,200);
y = linspace(-10,10,200);

    
[X,Y] = meshgrid(x,y);

% Define cylinder
a =1; %radius
phi = linspace(-pi,pi, 100);

[x_cldr,y_cldr] = pol2cart(phi,a);


% params
f=100;
k=0;
p0=1;
M=1000;
coords = 'cartesian';
t=0;
cyl_ctr = [0,0];
[p_i,p_s] = compute_cylinder_field(X,Y,a,f,k,p0,M,coords,t, cyl_ctr);


p_tot = p_s + p_i;