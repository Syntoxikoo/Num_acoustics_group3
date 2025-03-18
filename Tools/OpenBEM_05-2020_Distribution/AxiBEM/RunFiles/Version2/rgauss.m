function [bp,wf]=rgauss(n_intervals,nperint)
% function [bp,wf]=rgauss(n_intervals,nperint)
%
% Returns the weights and base points for the repeated Gauss
% numerical integration formula
% This is ideal for rapidly oscillating integrands
% INPUT
%     n_intervals : Number of intervals
%     nperint     : Number of integration points in each interval

% By msj 000222

[bpg,wfg]=gaussrule(nperint);

bp=bpg*ones(1,n_intervals)./n_intervals+ ...
   ones(nperint,1)*((0:n_intervals-1)-(n_intervals-1)/2)*2/n_intervals;
wf=wfg*ones(1,n_intervals)./n_intervals;

[bpi,bpj,bp]=find(bp-3);
[wfi,wfj,wf]=find(wf-3);
bp=bp+3;
wf=wf+3;

if any(bpi~=wfi) | any(bpj~=wfj)
   [bpi wfi bpj wfj]
   error('I!=J');
end    

if abs(sum(wf)-2)>1e-4 | find(abs(bp)>=1)
   sum(wf)
   bp
   error('Error in integration formula');
end
