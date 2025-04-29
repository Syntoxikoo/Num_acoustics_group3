addpath(genpath("tools"))
addpath(genpath("2DBEM"))
clear;

% calculate ground truth:
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

fr = 8000;
% el_wl=10*fr/c;
% [fr,GT_F,rr,theta] = flushed2DBEM(el_wl,fr);
% [~,GT_UF,~,~] = unflushed2DBEM(el_wl,fr);
% save("data/2DBEM_GT.mat","theta","rr","GT_UF","fr","GT_F");
load("data/2DBEM_GT.mat");
%% calculate convergence 
el_wl_vec = 1:0.25:6;
el_wl_vec=el_wl_vec*fr/c;
error_F = zeros(1,length(el_wl_vec));
error_UF = zeros(1,length(el_wl_vec));

for intI = 1:length(el_wl_vec)
    [~,F,~,~] = flushed2DBEM(el_wl_vec(intI),fr);
    [~,UF,~,~] = unflushed2DBEM(el_wl_vec(intI),fr);

    error_F(intI) = sum(abs(GT_F(1:length(theta))-F(1:length(theta))).^2./abs(GT_F(1:length(theta))).^2,1);
    error_UF(intI) = sum(abs(GT_UF(1:length(theta)) - UF(1:length(theta),:)).^2./abs(GT_UF(1:length(theta))).^2,1);
end
close all;
%%

semilogy(el_wl_vec,error_F,'DisplayName','flanged');
hold on;
semilogy(el_wl_vec,error_UF,'DisplayName','unflanged');
hold off;
grid on;
legend;
saveas(gcf, "project/figures/bem_conv1.svg")