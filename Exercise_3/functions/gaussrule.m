function [bp,wf]=gaussrule(n)
    % n: order, bp: base points, wf: weight factors
    u=(1:n-1)./sqrt((2*(1:n-1)).^2-1);
    [vc,bp]=eig(diag(u,-1)+diag(u,1));
    [bp,k]=sort(diag(bp));
    wf=2*vc(1,k)'.^2;
end