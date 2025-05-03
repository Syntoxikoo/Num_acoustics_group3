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
fp_theta = 1:5:360;
fp_theta = fp_theta';
fr = [1000;2000;4000;8000];
el_wl = 6*max(fr)/c;

%% Calculate
[fr,p_fieldUF,rr,fp_thetaUF,xybUF,topologyUF] = unflushed2DBEM(el_wl,fr,fp_theta);
[~,p_fieldR,~,fp_thetaR,xybR,topologyR] = unflushedRounded2DBEM(el_wl,fr,fp_theta);
[~,p_fieldSl,~,fp_thetaSl,xybSl,topologySl] = slit2DBEM(el_wl,fr,fp_theta);

%% plot

setup_RESULTs;


