function I=simpson(f,a,b,n)
% compute an integral using the simpson rule
% f = a function (form y=@(x) sqrt(x))
% a = lower end
% b = upper end
% n = number of intervals should be even!


h=(b-a)/n;
I1=0;
I2=0;

for i=1:n/2-1
    I1=I1+f(a+2*i*h);
end

for i=1:n/2
    I2=I2+f(a+(2*i-1)*h);
end

I=h/3*(f(a)+2*I1+4*I2+f(b));

end