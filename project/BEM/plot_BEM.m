% plotting the to compare

load("data/BEM_f.mat");
load("data/BEM_uf.mat");

figure;

theta_full = [theta; pi + theta];

corder = colororder;

for i = 1:length(fr)
    pF = p_fieldF(:,i);
    spl_valuesF = 20*log10(abs(pF(1:length(rr)))/20e-6);
    normalized_splF = spl_valuesF - max(spl_valuesF);
    spl_fullF = [normalized_splF; flip(normalized_splF)];

    p = p_fieldUF(:,i);
    spl_valuesUF = 20*log10(abs(p(1:length(rr)))/20e-6);
    normalized_spl_UF = spl_valuesUF - max(spl_valuesUF);
    spl_fullUF = [normalized_spl_UF; flip(normalized_spl_UF)];

    plot_F = polarplot(theta_full, spl_fullF, 'LineWidth', 2, 'DisplayName', [num2str(fr(i)) ' Hz'],"Color",corder(i,:));
    hold on;
    plot_UF = polarplot(theta_full, spl_fullUF, 'LineStyle','--','LineWidth', 2, 'HandleVisibility', 'off',"Color",corder(i,:));
    
end

pax = gca;
pax.ThetaZeroLocation = "top";
pax.ThetaDir = "clockwise";

grid on;
title('Directivity Pattern Comparison');
rlim([-60, 0]);
rticks(-60:6:0);
thetaticks(-180:15:180);
thetalim([-180 180]);
lgd = legend;
lgd.Position = [0.85, 0.1, 0.1, 0.1];
hold off;
