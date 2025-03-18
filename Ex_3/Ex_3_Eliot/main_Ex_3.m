clear;
clc;
close all;
addpath(genpath('function'))
addpath(genpath('scripts'))
%% part 1 : verif Num Int


f = @(x) sqrt(x);
a = 1; b = 4;
N= 1000;

I_r = riemann_sum(f,a,b,N);
I_t = trapezoidal_sum(f,a,b,N);
I_s = simpson_sum(f,a,b,N);
I_g = gauss_sum(f,a,b,N);

disp("Integral with Riemann sum : I="+I_r)
disp("Integral with trapezoid sum : I="+I_t)
disp("Integral with simpson sum : I="+I_s)
disp("Integral with Gauss sum : I="+I_g)

%% part 2 : Evaluate convergence
f = @(x) sqrt(x);
a = 1; b = 4;
N = 2:4:1000;
analytic = 14/3;
err_r = zeros(length(N),1);
err_t = zeros(length(N),1);
err_s = zeros(length(N),1);
err_g = zeros(length(N),1);
for i = 1:length(N)
    disp(i)
    I_r = riemann_sum(f,a,b,N(i));
    I_t = trapezoidal_sum(f,a,b,N(i));
    I_s = simpson_sum(f,a,b,N(i));
    I_g = gauss_sum(f,a,b,N(i));
    err_r(i) = abs(analytic - I_r);
    err_t(i) = abs(analytic - I_t);
    err_s(i) = abs(analytic - I_s);
    err_g(i) = abs(analytic - I_g);
end


figure(1)
loglog(N, err_r, 'LineWidth',1); hold on; grid on
loglog(N, err_t, 'LineWidth',1);
loglog(N, err_s, 'LineWidth',1);
loglog(N, err_g, 'LineWidth',1);
hold off;
xlim([2 1000])
legend('\epsilon riemann','\epsilon trapezoid','\epsilon simpson','\epsilon gauss')
xlabel('number of evaluation N')
ylabel("absolute error")
title("Convergence plot for differents numerical integration method : lb = 1 up= 4")

%% part 3 : for another int band

f = @(x) sqrt(x);
a = 0; b = 1;
N = 2:4:1000;
analytic = 2/3;
err_r = zeros(length(N),1);
err_t = zeros(length(N),1);
err_s = zeros(length(N),1);
err_g = zeros(length(N),1);
for i = 1:length(N)
    disp(i)
    I_r = riemann_sum(f,a,b,N(i));
    I_t = trapezoidal_sum(f,a,b,N(i));
    I_s = simpson_sum(f,a,b,N(i));
    I_g = gauss_sum(f,a,b,N(i));
    err_r(i) = abs(analytic - I_r);
    err_t(i) = abs(analytic - I_t);
    err_s(i) = abs(analytic - I_s);
    err_g(i) = abs(analytic - I_g);
end


figure(2)
loglog(N, err_r, 'LineWidth',1); hold on; grid on
loglog(N, err_t, 'LineWidth',1);
loglog(N, err_s, 'LineWidth',1);
loglog(N, err_g, 'LineWidth',1);
hold off;
xlim([2 1000])
legend('\epsilon riemann','\epsilon trapezoid','\epsilon simpson','\epsilon gauss')
xlabel('number of evaluation N')
ylabel("absolute error")
title("Convergence plot for differents numerical integration method\n lb = 0 up= 1")
