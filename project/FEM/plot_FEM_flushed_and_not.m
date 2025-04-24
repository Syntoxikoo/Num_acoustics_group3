
load("Result_FEM_unflushed.mat")
figure;
corder = colororder;
for ii =1: length(fr)

    Node_idx = cell2mat(N_arr(ii));

    theta = cell2mat(the_arr(ii));
    sortIdx = cell2mat(sor_arr(ii));
    p = cell2mat(p_arr(ii));
    spl_values = 20*log10(abs(p(Node_idx(sortIdx)))/20e-6);
    Nspl = spl_values - max(spl_values);
    polarplot(theta-pi/2, Nspl, 'LineWidth', 2, 'DisplayName', [num2str(fr(ii)) ' Hz'],Color=corder(ii,:));

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
    polarplot(theta-pi/2, Nspl, 'LineWidth', 2,"LineStyle", "--",'DisplayName', [num2str(fr(ii)) ' Hz'],Color=corder(ii,:));

    hold on;
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
legend('Location', 'best',NumColumns=2);
hold off;