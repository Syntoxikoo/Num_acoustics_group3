function I=gauss(f,a,b,n)
% compute an integral using the gaus rule
% f = a function (form y=@(x) sqrt(x))
% a = lower end
% b = upper end
% n = number of intervals should be even!

[bp,wf]=gaussrule(n);

bp = bp*(b-a)/2 + (b+a)/2;

I=0;

I=sum(wf.*f(bp))*(b-1)/2;

end