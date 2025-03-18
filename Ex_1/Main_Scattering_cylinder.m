clear all
close all
clc
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
%% Compute pressure around cylinder in far field
calc_press_polar
%% plot pressure around cylinder in far field
plot_press_polar
%% Anim pressure cart field
calc_press_anim