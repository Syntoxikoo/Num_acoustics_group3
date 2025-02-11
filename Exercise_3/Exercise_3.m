%% Exrcise 3 1)
clear;

a=1;
b=4;
n=100;

f= @(x) sqrt(x);

ILR=leftRiemann(f,a,b,n);
Itrap=trapez(f,a,b,n);
Isimp=simpson(f,a,b,n);
Igauss=gauss(f,a,b,n);

%% 2)
clear;

a=1;
b=4;
f= @(x) sqrt(x);

RES=14/3;

n=10:1:1000;

for i=1:1:length(n)
    ILR(i)=abs(RES-leftRiemann(f,a,b,i))/abs(RES);
    Itrap(i)=abs(RES-trapez(f,a,b,i))/abs(RES);
    Isimp(i)=abs(RES-simpson(f,a,b,i))/abs(RES);
    Igauss(i)=abs(RES-gauss(f,a,b,i))/abs(RES);
end

figure(2);

loglog(n,ILR);
hold on;
loglog(n,Itrap);
loglog(n,Isimp);
loglog(n,Igauss);
grid on;

legend('ILR','Itrap','Isimp','Igauss');

%% 3)
clear;

a=0;
b=1;
f= @(x) sqrt(x);

RES=14/3;

n=10:2:1000;

for i=1:1:length(n)
    ILR(i)=abs(RES-leftRiemann(f,a,b,i))/abs(RES);
    Itrap(i)=abs(RES-trapez(f,a,b,i))/abs(RES);
    Isimp(i)=abs(RES-simpson(f,a,b,i))/abs(RES);
    Igauss(i)=abs(RES-gauss(f,a,b,i))/abs(RES);
end

figure(2);

loglog(n,ILR);
hold on;
loglog(n,Itrap);
loglog(n,Isimp);
loglog(n,Igauss);
grid on;

legend('ILR','Itrap','Isimp','Igauss');