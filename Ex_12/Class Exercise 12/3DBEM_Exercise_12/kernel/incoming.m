function pI=incoming(srs,xyzb,k,varargin)

% pI=incoming(sources,xyzb,k)
%
% Computes the incoming field on the nodes.
% 
%  Input parameters (in SI units):
%  
%    -sources    : Sources information matrix. Each row is a source.
%                   Column 1:   if NaN, plane wave, else this is point source's
%                               xs coordinate.
%                   Column 2,3: theta (0 <= angle <= PI) and phi (0 <= angle <= 2*PI)
%                               for the plane wave, or point source's ys, zs coordinates.
%                   Column 4:   wave amplitude.
%                   Column 5:   initial phase (point source) or
%                               phase at origin (plane wave).
%    -xyzb       : Node coordinates. Matrix with M rows, x y z coordinates
%                  and the number of the body it belongs to.
%    -k          : Wavenumber (k=2*pi*f/c)
%
%  Output parameters:
%  
%    -pI         : Incoming field on the nodes (complex vector).
%                  A vector of zeros with lenght the number of nodes
%                  if the "sources" variable is empty.

% Vicente Cutanda 7-2001. 
% Correction VCH 2-2012.


if nargin>3
   disp('Too many input parameters in the ''incoming'' function. Extra input removed.')
end
   
pI=zeros(size(xyzb,1),1);
if ~isempty(srs)
   for ss=1:size(srs,1)
      if isnan(srs(ss,1)) % plane wave
         A=sin(srs(ss,2))*cos(srs(ss,3));
         B=sin(srs(ss,2))*sin(srs(ss,3));
         C=cos(srs(ss,2));
         % distance from the origin to the plane A(x-xi) + B(y-yi) + C(z-zi) = 0
         dps=(A*xyzb(:,1)+B*xyzb(:,2)+C*xyzb(:,3))./sqrt(A^2+B^2+C^2);
         pI=pI+srs(ss,4).*exp(-i*(k*dps+srs(ss,5)));
      else % point source
         % distance point source to nodes
         dps=sqrt((srs(ss,1)-xyzb(:,1)).^2+(srs(ss,2)-xyzb(:,2)).^2+(srs(ss,3)-xyzb(:,3)).^2);
         pI=pI+srs(ss,4)./dps.*exp(-i*(k*dps+srs(ss,5)));
      end
   end
end
