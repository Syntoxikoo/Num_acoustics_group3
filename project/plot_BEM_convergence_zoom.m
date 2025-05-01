% load data
clear;
load("2D_BEM/data/Conv_try12_zoomed_2.mat");
addpath plot_tools/

% calculate x axis vector
x_vec=el_wl_vec*c/fr;

nrows = 1;
ncols = 1;
heightScale = 0.7; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

tile1 = nexttile;

Leg(1) = semilogy(x_vec, error_UF(1,:), 'LineWidth', 1.5,Color="#181748"); hold on ; grid on 
%Leg(1) = semilogy(x_vec, error_UF(1,:), 'LineWidth', 1.5,Color="#810100"); grid on;

xlabel("Mesh density in elements per wavelength")
ylabel("mean error (pa)")
ylim([0 1])
leg1 = legend(Leg, "Recessed",NumColumns=2);
leg1.Layout.Tile = 'north';
title("Mesh convergence for BEM")
hold off
saveas(gcf,'figures/BEM_convergence_zoom.svg')