function [bp,wf]=nsingrule(n,xising,ndist)

% Returns the weights and base points for integration with
% a near singularity at 'xising'. Subinterval are integrated
% with the Gauss formula with 'n' points.
% 'ndist' is the normalised distance of node to element
%
% Details on this technique and its use with thin bodies
% and narrow gaps are given in:
% "Numerical transducer modelling", by Vicente Cutanda
% Proceedings of the 6th International Congress on Sound and Vibration
% Copenhaguen, July 1999


% By msj 990119
% Based on 'nsingint.m' by Vicente Cutanda
% N.B. Can be optimised like gaussrule

base=2;

if xising<-1 | xising>1
   error('Invalid xising value in local function nsingrule');
end

divs1=[1+xising];
while divs1(length(divs1))>ndist
   divs1=[divs1 divs1(length(divs1))/base];
end

divs2=[1-xising];
while divs2(1)>ndist
   divs2=[divs2(1)/base divs2];
end

if abs(xising)>=1
   divs=[xising-divs1 xising+divs2];
else
   divs=[xising-divs1 xising xising+divs2];
end

bp=[]; wf=[];
for idiv=2:length(divs)
   ldiv=abs(divs(idiv)-divs(idiv-1));
   [bpg,wfg]=gaussrule(n);
   bp=[bp; bpg*ldiv/2+(divs(idiv)+divs(idiv-1))/2];
   wf=[wf; wfg*ldiv/2];
end

if any(abs(bp)>=1) error('abs(bp) > 1'); end
if abs(sum(wf)-2)>eps*1e3 error(sprintf('sum(wf) = %g',sum(wf))); end
