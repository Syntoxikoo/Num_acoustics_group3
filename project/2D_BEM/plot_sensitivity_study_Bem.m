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

Leg(1) = semilogx(fr, 20 * log10(abs(pF_1)/2e-5), 'LineWidth', 1.5,Color="#181748"); hold on ; grid on 
Leg(2) = semilogx(fr, 20 * log10(abs(pF_2)/2e-5), 'LineWidth', 1.5,Color="#810100");

xlabel("Frequency (Hz)")
ylabel("Level (dB re 20 \mu Pa)")
ylim([40 80])
leg1 = legend(Leg, "10 cm Recessed","15 cm Recessed",NumColumns=2);
leg1.Layout.Tile = 'north';
title("Effect of change in piston depth " +rad2deg(theta)+"° from the axis")
hold off
% saveas(gcf,'figures/sensitivity_study_BEM_recessing.svg')