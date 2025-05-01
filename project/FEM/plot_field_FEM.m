load("Result_FEM_field.mat")
load("cmap.mat")
addpath("/Applications/MATLAB_R2024a.app/toolbox/pde")
phi = linspace(0,2*pi,50);

%% Fig 1

figure;
p = cell2mat(p_arr(1));
pdeplot(model1, 'XYData', abs(p), 'Contour', 'off', 'ColorMap', smoothColormap);
% clim([-0.8 0.8])

% xlim([-1.2 1.2]);
% ylim([-1.2 1.2]);


title("flushed")

saveas(gcf,"project/figures/FEM_field_flushed2k.svg")

figure;
set(gcf, 'Color', 'white');
for i = 1:length(phi)
    pdeplot(model1, "XYData", real(p*exp(1j*phi(i))),'ColorMap', smoothColormap);
    clim([-0.8 0.8])
    drawnow; % Force plot update

    pause(0.1); % Adjust speed
    frame(i) = getframe(gcf); % Capture frames for video export
end


v = VideoWriter('project/figures/FEM_anim_flushed.mp4', 'MPEG-4');
open(v);
writeVideo(v, frame);
close(v);

figure;
p = cell2mat(p_arr(2));
pdeplot(model2, 'XYData', abs(p), 'Contour', 'off', 'ColorMap', smoothColormap)
% clim([-0.8 0.8])


% xlim([-1.2 1.2]);
% ylim([-1.2 1.2]);


title("Recessed")
saveas(gcf,"project/figures/FEM_field_recessed2k.svg")


figure;
set(gcf, 'Color', 'white');
for i = 1:length(phi)
    pdeplot(model2, "XYData", real(p*exp(1j*phi(i))),'ColorMap', smoothColormap);
    clim([-0.8 0.8])
    drawnow; % Force plot update
    pause(0.1); % Adjust speed
    frame(i) = getframe(gcf); % Capture frames for video export
end


v = VideoWriter('project/figures/FEM_anim_recessed.mp4', 'MPEG-4');
open(v);
writeVideo(v, frame);
close(v);




% plot(xyb_chief(:,1),xyb_chief(:,2),"Color","black","Marker","square","Linestyle","none","MarkerSize", 5);
% xlabel(tiled, "x (m)")
% ylabel(tiled, "y (m)")
% title( M_high + " nodes, "+ N_high +" elements, "+ Nfield_high + " Field p.")
% grid on;
% hold off;

% leg = legend(Leg, {"Mesh body", "Field point", "Velocity"}, 'NumColumns', 3); 
% leg.Layout.Tile = 'north'; 

% saveas(gcf,'project/figures/mesh_size.svg')



