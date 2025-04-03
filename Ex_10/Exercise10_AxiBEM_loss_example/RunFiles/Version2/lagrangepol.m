function p=lagrangepol(n)

% Lagrange polinoms, order n

x = linspace(-1,1,n+1);
p = [];
for ii = 1:length(x)
    pp = 1;
    for jj = 1:length(x)
        if jj~=ii
            pp = conv(pp, [1 -x(jj)])/(x(ii)-x(jj));
        end
    end
    p = [p; pp];
end