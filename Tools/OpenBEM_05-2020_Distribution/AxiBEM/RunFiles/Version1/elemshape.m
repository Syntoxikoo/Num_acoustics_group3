function [psi, rhoq, zq, nrho, nz, dpsids, d2psids2]=elemshape(elknrzb,x)

% Calculate coordinates and normalvector at 'x'
% in an isoparametric element. (either linear or quadratic).
%
%  Input:
%    x:       real vector containing the ordinates for the
%             integration points (IP's)
%    elknrzb: real matrix, one row for each node in the element
%             each row contains (rho,z,body) for the node
%
%  Output:
%    psi:     Real matrix containing one shape function in each
%             column (no. of columns=no. of nodes in an element)
%    rhoq:    Real column vector containing the rho-coordinates
%             of the element for each IP
%    zq:      Real column vector containing the z-coordinates of
%             the element for each IP
%    nrho,nz: The column vectors containing respectively the rho
%             and z-coordinates of the normal vector.
%             Note that the normal vector does not have unit length.
%             This is done to avoid division by the jacobian which
%             occasionally is zero. Subsequent multiplication must
%             therefore be omitted.
%    dpsids:  Real matrix containing the derivative of 'psi' wrt
%             the tangent
%    d2psids2:Real matrix containing the second derivative of 'psi'
%             wrt. the tangent

[nknel,dummy]=size(elknrzb);

p = lagrangepol(nknel-1);
%p = bezier(nknel-1);

for aa = 1:size(p,1)
    psi(:,aa) = polyval(p(aa,:), x);
    dN(:,aa) = polyval(polyder(p(aa,:)),x);
end

% if nknel==2
%    % Linear Shape functions:
%    psi(:,1)=0.5*(1-x);
%    psi(:,2)=0.5*(1+x);
% 
%    % Linear Shape function derivatives
%    dN(:,1)=-0.5*ones(size(x));
%    dN(:,2)= 0.5*ones(size(x));
%    if nargout>6
%       d2N=zeros(size(dN));
%    end
% elseif nknel==3
%    % Quadratic Shape functions:
%    psi(:,1)=0.5*x.*(x-1);
%    psi(:,2)=1-x.^2;
%    psi(:,3)=0.5*x.*(x+1);
% 
%    % Quadratic Shape function derivatives
%    dN(:,1)=x-0.5;
%    dN(:,2)=-2*x;
%    dN(:,3)=x+0.5;
%    if nargout>6
%       d2N=ones(size(dN));
%       d2N(:,2)=-2*d2N(:,2);
%    end
% elseif nknel==4
%    p4 = [0.5625 0.5625 -0.0625 -0.0625];
%    p3 = [-1.6875 -0.5625 1.6875 0.5625];
%    p2 = [1.6875 -0.5625 -1.6875 0.5625];
%    p1 = [-0.5625 0.5625 0.0625 -0.0625];
%    
%    psi(:,1) = polyval(p1,x);
%    psi(:,2) = polyval(p2,x);
%    psi(:,3) = polyval(p3,x);
%    psi(:,4) = polyval(p4,x);
%     
%    dN(:,1) = polyval(polyder(p1),x);
%    dN(:,2) = polyval(polyder(p2),x);
%    dN(:,3) = polyval(polyder(p3),x);
%    dN(:,4) = polyval(polyder(p4),x);
%    
% elseif nknel==5   
%    p5 = [0.6667 0.6667 -0.1667 -0.1667 0];
%    p4 = [-2.6667 -1.3333 2.6667 1.3333 0];
%    p3 = [4 0 -5 0 1];
%    p2 = [-2.6667 1.3333 2.6667 -1.3333 0];
%    p1 = [0.6667 -0.6667 -0.1667 0.1667 0];
%     
%    psi(:,1) = polyval(p1,x);
%    psi(:,2) = polyval(p2,x);
%    psi(:,3) = polyval(p3,x);
%    psi(:,4) = polyval(p4,x);
%    psi(:,5) = polyval(p5,x);
%     
%    dN(:,1) = polyval(polyder(p1),x);
%    dN(:,2) = polyval(polyder(p2),x);
%    dN(:,3) = polyval(polyder(p3),x);
%    dN(:,4) = polyval(polyder(p4),x);
%    dN(:,5) = polyval(polyder(p5),x);
%     
% else
%    error('Only linear & quadratic & cubic & forth order? shape functions are implemented');
%    % Compute general shape function
%    % Still needs to be vectorised
%    [nip,dummy]=size(x);
%    psi=ones(nip+1,nknel);
%    for i=1:nknel
%       xnode(i)=(i-(nknel+1)/2)*2/(nknel-1+eps^2);
%    end
%    for i=1:nknel
%       x(nip+1)=xnode(i);
%       for k=1:nknel
%          if k~=i
%             psi(:,i)=psi(:,i)*(x-xnode(k));
%          end
%       end
%       psi(1:nip,i)=psi(1:nip,i)/psi(nip+1,i);
%    end
%    psi=psi(1:nip,:);
%    
%    psi
%    pause
%    % dN is still missing
% end

rhoq=psi*elknrzb(:,1);
zq=psi*elknrzb(:,2);

drho=dN*elknrzb(:,1);
dz=dN*elknrzb(:,2);

% The normal vector points into the domain, therefore the sign
% is changed if the domain is interior to the body.(VC)
nrho=-dz*sign(elknrzb(1,3));
nz=drho*sign(elknrzb(1,3));

if nargout>5
   dpsids=dN./(sqrt(drho.^2+dz.^2)*ones(1,nknel));
   if nargout>6
      d2rho=d2N*elknrzb(:,1);
      d2z=d2N*elknrzb(:,2);
      d2psids2=(d2N-dN.*(((drho.*d2rho+dz.*d2z)./(drho.^2+dz.^2))*ones(1,nknel)))./((drho.^2+dz.^2)*ones(1,nknel));
   end
end
