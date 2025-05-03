% clear

% nrows = 1;
% ncols = 2;
% heightScale = 0.7; % Adjust height scaling if needed
% Nleg = 1;

% [columnwidth, ~] = get_widths();
% height = get_height() * heightScale; 
% fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
% tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
% corder = ["#181748","#810100","#575579","#7F7F7F"];


% tile1 = nexttile;


% load("Result_FEM_unflushed.mat")
% Node_idx = cell2mat(N_arr(2));
% theta = cell2mat(the_arr(2));
% sortIdx = cell2mat(sor_arr(2));
% p = cell2mat(p_arr(2));
% spl_values = 20*log10(abs(p(Node_idx(sortIdx)))/20e-6);
% Nspl = spl_values - max(spl_values);
% Leg(1) = polarplot(theta-pi/2, Nspl, 'LineWidth', 1.5,Color=corder(1));

% hold on;

% load("2DBEM_uf.mat")
% p = p_fieldUF(:,2);
% spl_values = 20*log10(abs(p(1:length(rr)))/20e-6);
% Nspl = spl_values - max(spl_values);
% Leg(2) = polarplot(theta, Nspl, 'LineWidth', 1.5,"Linestyle",'--',Color=corder(2));


% pax = gca;
% pax.ThetaZeroLocation = "top";
% pax.ThetaDir = "clockwise";
% title( "2D FEM")
% grid on;
% rlim([-60, 0]);
% rticks(-60:6:0);
% thetaticks(-180:15:180);
% thetalim([-180 180]);
% leg1 = legend(Leg,'Location', 'best',NumColumns=4);
% leg1.Layout.Tile = 'north';
% hold off;

% tile2 = nexttile;

% load("2DBEM_uf.mat")

% load("Result_FEM_unflushed.mat")
% Node_idx = cell2mat(N_arr(2));
% theta1 = cell2mat(the_arr(2));
% sortIdx = cell2mat(sor_arr(2));
% p1 = cell2mat(p_arr(2));
% p1 =p1(Node_idx(sortIdx));
% load("2DBEM_uf.mat")
% theta2 = theta(1:end-1);
% spl_values = 20*log10(abs(p1)/20e-6);
% Nspl1 = spl_values - max(spl_values);
% p = p_fieldUF(:,2);
% spl_values2 = 20*log10(abs(p(1:length(theta2)))/20e-6);
% Nspl2 = spl_values2 - max(spl_values2);
% Nspl = interp1(theta1,Nspl1,theta2);
% Nspl = circshift(Nspl,-90);
% % plot(rad2deg(theta2), Nspl - Nspl2, 'LineWidth', 1.5,Color=corder(3));
% polarplot(theta2, Nspl, 'LineWidth', 1.5,"Linestyle",'-',Color=corder(1)); hold on
% polarplot(theta2, Nspl2, 'LineWidth', 1.5,"Linestyle",'--',Color=corder(2));
% pax = gca;
% pax.ThetaZeroLocation = "top";
% pax.ThetaDir = "clockwise";
% title( "2D FEM")
% grid on;
% rlim([-60, 0]);
% rticks(-60:6:0);
% thetaticks(-180:15:180);
% thetalim([-180 180]);
% leg1 = legend(Leg,'Location', 'best',NumColumns=4);
% leg1.Layout.Tile = 'north';
% hold off;

% title( "2D BEM")

% % ----------------


% hold off;

% % saveas(gcf,'project/figures/2DBEM_FEM_comp.svg')
clear

nrows = 1;
ncols = 2;
heightScale = 0.7; % Adjust height scaling if needed

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = ["#181748","#810100","#575579","#7F7F7F"];

% Load data once to avoid redundant loading
load("Result_FEM_unflushed.mat")
load("2DBEM_uf.mat")

% Tile 1: Compare directivity patterns
tile1 = nexttile;

% FEM data
Node_idx = cell2mat(N_arr(2));
theta_fem = cell2mat(the_arr(2));
sortIdx = cell2mat(sor_arr(2));
p_fem = cell2mat(p_arr(2));
p_fem = p_fem(Node_idx(sortIdx));
spl_fem = 20*log10(abs(p_fem)/20e-6);
Nspl_fem = spl_fem - max(spl_fem);

% BEM data
theta_bem = theta(1:end-1);  % Using all data points except the last one
p_bem = p_fieldUF(:,2);
p_bem = p_bem(1:length(theta_bem));  % Match length with theta_bem
spl_bem = 20*log10(abs(p_bem)/20e-6);
Nspl_bem = spl_bem - max(spl_bem);

% Plot both directivity patterns
Leg(1) = polarplot(theta_fem-pi/2, Nspl_fem, 'LineWidth', 1.5, 'Color', corder(1));
hold on;
Leg(2) = polarplot(theta_bem, Nspl_bem, 'LineWidth', 1.5, 'Linestyle', '--', 'Color', corder(2));

pax = gca;
pax.ThetaZeroLocation = "top";
pax.ThetaDir = "clockwise";
title("Directivity Patterns Comparison")
grid on;
rlim([-60, 0]);
rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);
leg1 = legend(Leg, {'FEM', 'BEM'}, 'Location', "best", 'NumColumns', 1);
hold off;

% Tile 2: Show difference between directivity patterns
tile2 = nexttile;

% Interpolate FEM data to BEM angles for comparison
Nspl_fem_interp = interp1(theta_fem, Nspl_fem, theta_bem);
% Align patterns by shifting as needed
Nspl_fem_interp = circshift(Nspl_fem_interp, -90);

% Calculate the difference
pattern_difference = Nspl_fem_interp - Nspl_bem;

% Plot the difference
plot(rad2deg(theta_bem)-180, pattern_difference, 'LineWidth', 2, 'Color',"#000");
xlabel('Angle (degrees)');
ylabel('Difference (dB)');
title("Difference (FEM - BEM)");
grid on;
xlim([-180 180]);
xticks(-180:30:180);
ylim([-1.5,1.5])

% Add a zero reference line
hold on;
plot([0 360], [0 0], 'k--', 'LineWidth', 0.5);
hold off;

% Adjust overall layout
sgtitle('Comparison of 2D FEM vs 2D BEM', 'FontWeight', 'bold');

% Uncomment to save
saveas(gcf,'project/figures/2DBEM_FEM_comp.svg')