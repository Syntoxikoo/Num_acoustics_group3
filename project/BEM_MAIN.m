addpath(genpath("Axi_BEM"));
addpath(genpath("data"));
addpath(genpath("RunFiles"));

close all;
clear;

%% calculate Box
flushed_piston_AxiBEM;
clear;

%%
unflushed_piston_AxiBEM;
clear;

%% plot 
plot_AxiBEM;