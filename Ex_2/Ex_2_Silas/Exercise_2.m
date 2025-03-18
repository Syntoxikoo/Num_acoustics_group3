% Exercise 2
clear;
% first one

x=linspace(-1e-6,1e-6,100);


x=single(x);

y=sqrt(x+1)-1;

yfix=x./(sqrt(x+1)+1);

figure(1);
subplot(1,3,1);
plot(x,y);
hold on;
plot(x,yfix,"LineStyle","--");
hold off;
grid on;

% second one
clear;
x=linspace(-3,3,1000);
y=x+1e-6;

x=single(x);
y=single(y);

z=sin(x)-sin(y);

zfix=2.*cos((x+y)/2).*sin((x-y)/2);

subplot(1,3,2);
plot(x,z);
hold on;
plot(x,zfix,"LineStyle","--");
hold off;
grid on;

% last one
clear;
x=linspace(0,2,1000);

y=x-1e-6;

x=single(x);
y=single(y);
z2=x.^2-y.^2;

z2fix=(x+y).*(x-y);

grid on;

subplot(1,3,3);
plot(x, z2);
hold on;
plot(x,z2fix,"LineStyle","--");
hold off;
grid on;

