function I=leftRiemann(f,a,b,n)
% compute an integral using the left riemann rule
% f = a function (form y=@(x) sqrt(x))
% a = lower end
% b = upper end
% n = number of intervals

h=(b-a)/n;
I=0; 

for i=0:n-1
    I=I+f(a+i*h);
end

I=h*I;

end