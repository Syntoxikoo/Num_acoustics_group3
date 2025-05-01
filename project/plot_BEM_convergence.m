% load data
clear;
load("2D_BEM/data/Conv_try14_eliot_zoom.mat");
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

Leg(1) = semilogy(x_vec, error_F(1,:), 'LineWidth', 1.5,Color="#181748"); hold on ; grid on 
Leg(2) = semilogy(x_vec, error_UF(1,:), 'LineWidth', 1.5,Color="#810100");

xlabel("Mesh density in elements per wavelength")
ylabel("mean error (pa)")
ylim([1e-2 12])
xlim([1 20])
leg1 = legend(Leg,"flushed", "Recessed",NumColumns=2);
leg1.Layout.Tile = 'north';
title("Mesh convergence for BEM")
hold off
saveas(gcf,'figures/BEM_convergence.svg')

