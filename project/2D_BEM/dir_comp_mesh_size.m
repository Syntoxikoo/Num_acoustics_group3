
nrows = 1;
ncols = 2;
heightScale = 0.7; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

% title( "Directivity of the Loudspeaker for 2 spacial resolution")
tile1 = nexttile;


spl_values = 20*log10(abs(pF_low(1:length(rr_low)))/20e-6);
Nspl_low = spl_values - max(spl_values);
Leg(1)=polarplot(theta_low, Nspl_low, 'LineWidth', 1.5,Color="#181748");
hold on;

spl_values = 20*log10(abs(pF_low2(1:length(rr_low2)))/20e-6);
Nspl_low2 = spl_values - max(spl_values);
Leg(2)=polarplot(theta_low2, Nspl_low2, 'LineWidth', 1.5,Color="#810100");

pax = gca;
pax.ThetaZeroLocation = "top";
pax.ThetaDir = "clockwise";
title( "1 Elements per \lambda")
grid on;
rlim([-60, 0]);
rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);

hold off;

tile2 = nexttile;

spl_values = 20*log10(abs(pF_high(1:length(rr_high)))/20e-6);
Nspl_high = spl_values - max(spl_values);
polarplot(theta_high, Nspl_high, 'LineWidth', 1.5,Color="#181748");
hold on;
spl_values = 20*log10(abs(pF_high2(1:length(rr_high2)))/20e-6);
Nspl_high2 = spl_values - max(spl_values);
polarplot(theta_high2, Nspl_high2, 'LineWidth', 1.5,Color="#810100");



pax = gca;
pax.ThetaZeroLocation = "top";
pax.ThetaDir = "clockwise";
title("6 Elements per \lambda")
grid on;
rlim([-60, 0]);
rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);

leg1 = legend(Leg, "f= 2000 Hz","f= 8000 Hz",NumColumns=2);
leg1.Layout.Tile = 'south';
hold off;

saveas(gcf,'project/figures/dir_comp_mesh_size.svg')