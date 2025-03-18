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

t= 0;
cyl_ctr = [0,0];
[p_i,p_s] = compute_cylinder_field(X,Y,a,f,k,p0,M,coords,t, cyl_ctr);

t=linspace(0,1,100);
p_tot = p_s + p_i;
if length(t)>1
    p_tmp = zeros(size(p_tot,1),size(p_tot,2),length(t));
    for ii = 1: length(t)
        p_tmp(:,:,ii) = p_tot.* exp(-1j*2*pi*f*t(ii));
    end
    p_tot = p_tmp;
end
r= sqrt(X.^2+Y.^2);
mask = find(r <=a);
p_tot(mask) = NaN;
% 
hh = figure;
jj = 1;
xy=[X(1:end)' Y(1:end)'];
[Xc, Yc, Zc] = cylinder(a, 100);
Zc = Zc * max(max(abs(p_tot(:,:,1))));
Zc(Zc==0) = min(min(abs(p_tot(:,:,1))));

videoName = 'figures/anim_wave.mov';
vw = VideoWriter(videoName, 'MPEG-4');
vw.FrameRate = 30;  % 30 fps for smooth playback
open(vw);

startTime = tic;
elapsedTime = 0;

while elapsedTime <10
    figure(hh)
    surf(X,Y,real(p_tot(:,:,jj))); hold on
    xlabel('x, m'); ylabel('y, m');zlabel('p, Pa');
    axis([min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2)) min(real(reshape(p_tot,[],1))) max(real(reshape(p_tot,[],1)))])
    view(-30,40)
    surf(Xc, Yc, Zc) ; hold off
    % drawnow;
    frame = getframe(hh);  
    writeVideo(vw, frame);
    if jj<100
        jj = jj+1;
    else
        jj=1;
    end
    elapsedTime = toc(startTime);
    disp(elapsedTime)
end