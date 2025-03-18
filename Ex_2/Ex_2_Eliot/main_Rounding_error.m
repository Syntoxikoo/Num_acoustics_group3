clear variables;
close all;
clc;
addpath(genpath('scripts'));
addpath(genpath('functions'))

%% Rounding error 1
x = single(linspace(-1e-6, 1e-6, 1000));
y = sqrt(x+1) - 1;
y2= x./(sqrt(x+1)+1);
figure(1);
plot(x, y, '-', 'LineWidth',1);
hold on;
plot(x, y2, '-', 'LineWidth',1);
xlabel('x \approx 0');
ylabel('y = sqrt(x+1) - 1');
legend('formulation', 'reformulation');
title('Rounding error');
grid on;
%% Rounding error 2
x = single(linspace(-pi, pi, 10000));
y = single(linspace(-pi, pi, 10000))+1e-6;
z = sin(x) - sin(y);
z2 = sin((x-y)/2) .*2.*cos((x+y)/2); 
figure(2);
plot(x, z, '-', 'LineWidth',1); hold on;
plot(x, z2, '-', 'LineWidth',1);
grid on;
legend('formulation', 'reformulation');
xlabel('x \approx y');
ylabel('z = sin(x) - sin(y)');
title('Rounding error');
%% Rounding error 3
x = single(linspace(0, 2, 10000))+1e-6;
y = single(linspace(0, 2, 10000));
z = x.^2 - y.^2;
z2 = (x+y).*(x-y);
figure(3);
plot(x, z, '-', 'LineWidth',1); hold on;
plot(x, z2, '-', 'LineWidth',1);
xlabel('x \approx y');
ylabel('z = x^2 - y^2');
legend('formulation', 'reformulation');
grid on;
xlim([0 2]);
title('Rounding error');
