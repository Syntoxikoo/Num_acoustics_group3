function I = gauss_sum(f,a,b,N)
    % numericall integration using riemann sum method
    % Input: 
    %   - a : lower bound
    %   - b : upper bound
    %   - f : integrand
    %   - dx : step_size
    %   - method : 'left', 'right', 'mid'
    arguments
        f
        a = -1
        b = 1
        N = 1000
    end
    [bp,wf] = gaussRule(N);

    if a ~= -1 || b ~= 1
        bp = bp*(b-a)/2 + (b+a) /2;
        I = sum(wf .* f(bp));
        I = I* (b-a)/2;
    else
        I = sum(wf .* f(bp));
    end

end

function [bp,wf] = gaussRule(N)
    % Input :
    %   - N: order
    % Output : 
    %   - bp : base points
    %   - wf : weight factor
    u = (1:N-1)./sqrt((2*(1:N-1)).^2-1);
    [vc,bp] = eig(diag(u,-1)+diag(u,1));
    [bp,k] = sort(diag(bp));
    wf = 2 * vc(1,k)'.^2;
end