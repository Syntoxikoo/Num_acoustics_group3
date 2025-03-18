function I = simpson_sum(f,a,b,N,method)
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
        method = "1/3"
    end
    dx = (b-a)/N;
    
    if lower(method) == "1/3"
        I = dx/3 * (f(a) + f(b) + 2 * sum(f(a + 2*dx*(1:(N/2-1))) )+ 4 * sum(f(a + dx*(2*(1:(N/2))-1)) ));
    end
end
