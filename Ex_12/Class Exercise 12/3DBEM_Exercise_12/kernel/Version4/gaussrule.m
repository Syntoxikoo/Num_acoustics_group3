function [bp,wf]=gaussrule(n)

% Returns the weights and base points for the Gauss numerical integration
% formula with 'n' points

% Gauss base points and weight factors calculation taken from
% a user posted function in the Mathworks website:
% Concepts on page 93 of
% 'Methods of Numerical Integration' by
% Philip Davis and Philip Rabinowitz yield
% the base points and weight factors.
%
%          Howard Wilson
%          Department of Engineering Mechanics
%          University of Alabama
%          Box 870278
%          Tuscaloosa, Alabama 35487-0278
%          Phone 205 348-1617
%          Email address: HWILSON @ UA1VM.UA.EDU


persistent GAUSSRULE_BASE_POINTS GAUSSRULE_WEIGHTS

if n<1
   bp=[];
   wf=[];
elseif n<=size(GAUSSRULE_WEIGHTS) & GAUSSRULE_WEIGHTS(n,1)~=0
   bp=GAUSSRULE_BASE_POINTS(n,1:n)';
   wf=GAUSSRULE_WEIGHTS(n,1:n)';
else
   u=(1:n-1)./sqrt((2*(1:n-1)).^2-1);
   [vc,bp]=eig(diag(u,-1)+diag(u,1));
   [bp,k]=sort(diag(bp));
   wf=2*vc(1,k)'.^2;
   [oldn,dummy]=size(GAUSSRULE_WEIGHTS);
   if n>oldn
      GAUSSRULE_WEIGHTS=[GAUSSRULE_WEIGHTS zeros(oldn,n-oldn)];
      GAUSSRULE_BASE_POINTS=[GAUSSRULE_BASE_POINTS zeros(oldn,n-oldn)];
   end
   GAUSSRULE_WEIGHTS(n,1:n)=wf';
   GAUSSRULE_BASE_POINTS(n,1:n)=bp';
   sparse(GAUSSRULE_WEIGHTS);
   sparse(GAUSSRULE_BASE_POINTS);
end
