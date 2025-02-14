function I = trapezoidal_sum(f,a,b,N)
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
    end
    dx = (b-a)/N;
    x = (a:dx:b);
    
        % I = dx *((f(a)+f(b))/2 + sum(f(a+x(2:end-1))));
        I = dx * sum((f(x(1:end-1))+f(x(2:end)))/2);
end