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

 
u=(1:n-1)./sqrt((2*(1:n-1)).^2-1);
[vc,bp]=eig(diag(u,-1)+diag(u,1));
[bp,k]=sort(diag(bp));
wf=2*vc(1,k)'.^2; 


