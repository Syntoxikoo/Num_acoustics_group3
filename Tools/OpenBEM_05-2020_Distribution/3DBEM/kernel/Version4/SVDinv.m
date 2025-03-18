function [invMat,dS]=SVDinv(inMat,SVthld)

% [invMat,dS]=SVDinv(inMat,SVthld);
%
% Performs a matrix inversion of the matrix "inMat" using SVD. 
% It neglects any singular value smaller than "SVthld".
% dS: singular values

% Vicente Cutanda 08-2004

[U,S,V] = svd(inMat); % SVD for inverting Vy
dS=diag(S);
hh=find(dS>SVthld);
iS=zeros(length(dS),1);
iS(hh)=1./dS(hh);
invMat=V*diag(iS)*U';
