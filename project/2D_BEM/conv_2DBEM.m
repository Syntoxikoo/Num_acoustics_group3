addpath(genpath("tools"))
addpath(genpath("2DBEM"))
clear;
clc;

% calculate ground truth:
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

% define angles of observation
fp_theta = [33;55;78;135];
% fp_theta = fp_theta';
% calculate ground truth
fr = 4000;
el_wl=50*fr/c;
[fr,GT_F,rr,theta] = flushed2DBEM(el_wl,fr,fp_theta);
[~,GT_UF,~,~] = unflushed2DBEM(el_wl,fr,fp_theta);
%save("data/2DBEM_GT.mat","theta","rr","GT_UF","fr","GT_F");
% load("data/2DBEM_GT.mat");
%% calculate convergence 
el_wl_vec = 1:2:50;
el_wl_vec=el_wl_vec*fr/c;
error_F = zeros(2,length(el_wl_vec));
error_UF = zeros(2,length(el_wl_vec));
count = 1:1:length(el_wl_vec);

for intI = 1:length(el_wl_vec)

    % counter
    disp(['Simulation ' num2str(count(intI)) ' of ' num2str(count(end))]);

    [~,F,~,~] = flushed2DBEM(el_wl_vec(intI),fr,fp_theta);
    [~,UF,~,~] = unflushed2DBEM(el_wl_vec(intI),fr,fp_theta);
    close all;
    % Absolut
    error_F(1,intI) = mean(abs(GT_F-F)./abs(GT_F));
    error_UF(1,intI) = mean(abs(GT_UF-UF)./abs(GT_UF));
    % angle
    error_F(2,intI) = mean(angle(GT_F-F)./angle(GT_F));
    error_UF(2,intI) = mean(angle(GT_UF-UF)./angle(GT_UF));

    % error_F(intI) = sum(abs(GT_F-F).^2./abs(GT_F).^2,1);
    % error_UF(intI) = sum(abs(GT_UF - UF).^2./abs(GT_UF).^2,1);
    
end
close all;
%%

plot(el_wl_vec,error_F(1,:),'DisplayName','flanged');
hold on;
plot(el_wl_vec,error_UF(1,:),'DisplayName','unflanged');
hold off;
grid on;
legend;
%saveas(gcf, "project/figures/bem_conv1.svg")