addpath(genpath('scripts'))
addpath(genpath('functions'))
%% Compute pressure cartesian field
calc_press_2Dfield
%% Plot pressure 2D field
plot_press_2Dfield
%% Compute pressure along the cylinder surface
calc_press_cyl
%% plot pressure along the cylinder surface
plot_press_cyl_surf