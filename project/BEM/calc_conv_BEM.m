addpath(genpath("project"))
addpath(genpath("BEM"))
clear;

% calculate ground truth:
% Ambient conditions
pa = 101325;         % Static pressure (Pa)
t = 20;              % Temperature (ºC)
Hr = 50;             % Relative humidity (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

max_fr = 8000;
el_wl=6*max_fr/c;
% [fr,GT_F,rr,theta] = flushedBEM(el_wl);
% [~,GT_UF,~,~] = unflushedBEM(el_wl);
% save("data/BEM_GT.mat","theta","rr","GT_UF","fr","GT_F");
load("data/BEM_GT.mat");
%% calculate convergence 
el_wl_vec = 1:1:2;

error_F = zeros(length(el_wl_vec),4);
error_UF = zeros(length(el_wl_vec),4);

for intI = 1:length(el_wl_vec)
    [~,F,~,~] = flushedBEM(el_wl_vec(intI));
    [~,UF,~,~] = unflushedBEM(el_wl_vec(intI));

    error_F(intI,:) = sum(abs(GT_F(1:length(theta),:)-F(1:length(theta),:)).^2./abs(GT_F(1:length(theta),:)).^2,1);
    error_UF(intI,:) = sum(abs(GT_UF(1:length(theta),:) - UF(1:length(theta),:)).^2./abs(GT_UF(1:length(theta),:)).^2,1);
end

%%
semilogy(el_wl_vec,error_UF);