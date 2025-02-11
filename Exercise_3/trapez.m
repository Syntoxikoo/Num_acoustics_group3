function I=trapez(f,a,b,n)
% compute an integral using the trapez rule
% f = a function (form y=@(x) sqrt(x))
% a = lower end
% b = upper end
% n = number of intervals

h=(b-a)/n;
I=0; 

for i=1:n-1
    I=I+f(a+i*h);
end

I=h*(1/2*f(a)+I+1/2*f(b));
end