

nrows = 1;
ncols = 2;
heightScale = 0.7; % Adjust height scaling if needed
Nleg = 1;

[columnwidth, ~] = get_widths();
height = get_height() * heightScale; 
fig = figure("Position", [0, 0, columnwidth, height], "Units", "points");
tiled = tiledlayout(nrows, ncols, "TileSpacing", "tight", "Padding", "loose");
corder = colororder;

load("Result_FEM_unflushed.mat")

tile1 = nexttile;

for ii =1: length(fr)

    Node_idx = cell2mat(N_arr(ii));

    theta = cell2mat(the_arr(ii));
    sortIdx = cell2mat(sor_arr(ii));
    p = cell2mat(p_arr(ii));
    spl_values = 20*log10(abs(p(Node_idx(sortIdx)))/20e-6);
    Nspl = spl_values - max(spl_values);
    leg1(ii) = polarplot(theta-pi/2, Nspl, 'LineWidth', 1, 'DisplayName', [num2str(fr(ii)) ' Hz'],Color=corder(ii,:));

    hold on;
end

load("Result_FEM_flushed.mat")
for ii =1: length(fr)

    Node_idx = cell2mat(N_arr(ii));

    theta = cell2mat(the_arr(ii));
    sortIdx = cell2mat(sor_arr(ii));
    p = cell2mat(p_arr(ii));
    spl_values = 20*log10(abs(p(Node_idx(sortIdx)))/20e-6);
    Nspl = spl_values - max(spl_values);
    polarplot(theta-pi/2, Nspl, 'LineWidth', 1,"LineStyle", "--","Color",corder(ii,:),'HandleVisibility', 'off');

    hold on;
end


pax = gca;
pax.ThetaZeroLocation = "top";
pax.ThetaDir = "clockwise";
title( "2D FEM")
grid on;
rlim([-60, 0]);
rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);
leg1 = legend('Location', 'best',NumColumns=4);
leg1.Layout.Tile = 'south';
hold off;

tile2 = nexttile;

load("axi_BEM_Uf.mat")

for ii = 1:length(fr)
    theta = cell2mat(theta_arr(ii));
    theta_full = [theta; pi + theta];
    p = cell2mat(p_fieldF_arr(ii));
    spl_values = 20*log10(abs(p(1:length(cell2mat(theta_arr(ii)))))/20e-6);
    Nspl = spl_values - max(spl_values);
    spl_full = [Nspl; flip(Nspl)];
    Leg(1)=polarplot(theta_full, spl_full, 'LineWidth', 1, 'DisplayName', [num2str(fr(ii)) ' Hz'],Color=corder(ii,:));

    hold on;
end

load("axi_BEM_f.mat")
for ii = 1:length(fr)
    theta = cell2mat(theta_arr(ii));
    theta_full = [theta; pi + theta];
    p = cell2mat(p_fieldF_arr(ii));
    spl_values = 20*log10(abs(p(1:length(cell2mat(theta_arr(ii)))))/20e-6);
    Nspl = spl_values - max(spl_values);
    spl_full = [Nspl; flip(Nspl)];
    Leg(2) = polarplot(theta_full, spl_full, 'LineWidth', 1,"LineStyle", "--",'DisplayName', [num2str(fr(ii)) ' Hz'],Color=corder(ii,:));

    hold on;
end


pax = gca;
pax.ThetaZeroLocation = "top";
pax.ThetaDir = "clockwise";

grid on;
rlim([-60, 0]);
rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);

title( "Axisymetrical BEM")

% ----------------
h(1) = plot(NaN, NaN, 'k-', 'LineWidth', 1); 
h(2) = plot(NaN, NaN, 'k--', 'LineWidth', 1); 
leg = legend(h, {"unflushed piston", "flushed piston"}, 'NumColumns', 2); 
leg.Layout.Tile = 'north'; 

hold off;

saveas(gcf,'project/figures/axiBEM_FEM_comp.svg')