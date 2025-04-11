addpath(genpath("BEM"));
addpath(genpath("data"));
addpath(genpath("RunFiles"));

close all;
clear;

%% calculate Box
piston_baffle;
clear;

%%
unflushed_piston;
clear;

%% plot 
plot_BEM;