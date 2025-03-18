function I = riemann_sum(f,a,b,N,method)
    % numericall integration using riemann sum method
    % Input: 
    %   - a : lower bound
    %   - b : upper bound
    %   - f : integrand
    %   - dx : step_size
    %   - method : 'left', 'right', 'mid'
    arguments
        f
        a
        b
        N = 1000
        method = "left"
    end
    dx = (b-a)/N;
    x = (a:dx:b);
    if lower(method) == "left"
        I = sum(f(x)*dx);
    end
end