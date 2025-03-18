function [psi, xq, yq, zq, nx, ny, nz]=elemshape(elknxyzb,IP)

% [psi, xq, yq, zq, nx, ny, nz]=elemshape(elknxyzb,IP)
%
% Calculate shape functions, global coordinates and normal vector at
% a given set of local points IP=(e1i,e2i), where (e1,e2) are the local
% coordinates in an isoparametric element, either linear or quadratic.
%
%  Input:
%    IP:      real matrix containing the local coordinates for the 
%             integration points (IP's) e1 values in 1st column
%             e2 values in 2nd.
%    elknxyzb: real matrix, one row for each node in the element
%             each row contains (x,y,z,body) for the node.
%
%  Output:
%    psi:     Real matrix containing one shape function in each
%             column (no. of columns=no. of nodes in an element)
%    xq:      Real column vector containing the global x-coordinates
%             of the element for each IP.
%    yq:      Real column vector containing the global y-coordinates
%             of the element for each IP.
%    zq:      Real column vector containing the global z-coordinates
%             of the element for each IP.
%    nx,ny,nz: 
%             The column vectors containing respectively the global
%             x,y and z-coordinates of the normal vector.
%             Note that the normal vector does not have unit length.
%             This is done to avoid division by the jacobian which
%             occasionally is zero. Subsequent multiplication must
%             therefore be omitted.
%
% References: -Peter M. Juhl: "The Boundary Element Method for Sound
%             Field Calculations", Report No. 55, DTU 1993. Section 4.7.
%             -O. C. Zienkiewicz, R. L. Taylor: "The Finite Element Method"
%             4th Ed., Volume 1, Section 8.5.

% Peter M. Juhl 2000.

% Added more comments, local coordinate notation and faster formulation
% (Vicente Cutanda 4-2001).

% Light version with less calculation enabled, to be used in element
% subdivision for near-singular integrals. Vicente Cutanda Henriquez 12-2010.

nknel=size(elknxyzb,1); % number of nodes

switch nknel
    case {4}
        % Linear Shape functions:
        psi=[0.25*(1-IP(:,1)).*(1-IP(:,2)) ...
            0.25*(1+IP(:,1)).*(1-IP(:,2)) ...
            0.25*(1+IP(:,1)).*(1+IP(:,2)) ...
            0.25*(1-IP(:,1)).*(1+IP(:,2))];
        
        if nargout>4
            % Linear Shape function derivatives
            % Note that these elements are not necessary plane
            % (as opposed to those in my thesis). Hence, the
            % normal vector vary along the element
            dNde1=[-0.25*(1-IP(:,2)) ...
                0.25*(1-IP(:,2)) ...
                0.25*(1+IP(:,2)) ...
                -0.25*(1+IP(:,2))];
            
            dNde2=[-0.25*(1-IP(:,1)) ...
                -0.25*(1+IP(:,1)) ...
                0.25*(1+IP(:,1)) ...
                0.25*(1-IP(:,1))];
        end
        
    case {8}
        % Quadratic Shape functions.
        psi=[0.25*(1-IP(:,1)).*(1-IP(:,2)).*(-IP(:,1)-IP(:,2)-1) ...
            0.25*(1+IP(:,1)).*(1-IP(:,2)).*(IP(:,1)-IP(:,2)-1) ...
            0.25*(1+IP(:,1)).*(1+IP(:,2)).*(IP(:,1)+IP(:,2)-1) ...
            0.25*(1-IP(:,1)).*(1+IP(:,2)).*(-IP(:,1)+IP(:,2)-1) ...
            0.5*(1-IP(:,1).^2).*(1-IP(:,2)) ...
            0.5*(1-IP(:,2).^2).*(1+IP(:,1)) ...
            0.5*(1-IP(:,1).^2).*(1+IP(:,2)) ...
            0.5*(1-IP(:,2).^2).*(1-IP(:,1))];
        
        if nargout>4
            % Quadratic Shape functions derivatives.
            dNde1=[0.25*(1-IP(:,2)).*(2*IP(:,1)+IP(:,2)) ...
                0.25*(1-IP(:,2)).*(2*IP(:,1)-IP(:,2)) ...
                0.25*(1+IP(:,2)).*(2*IP(:,1)+IP(:,2)) ...
                0.25*(1+IP(:,2)).*(2*IP(:,1)-IP(:,2)) ...
                -IP(:,1).*(1-IP(:,2)) ...
                0.5*(1-IP(:,2).^2) ...
                -IP(:,1).*(1+IP(:,2)) ...
                -0.5*(1-IP(:,2).^2)];
            
            dNde2=[0.25*(1-IP(:,1)).*(2*IP(:,2)+IP(:,1)) ...
                0.25*(1+IP(:,1)).*(2*IP(:,2)-IP(:,1)) ...
                0.25*(1+IP(:,1)).*(2*IP(:,2)+IP(:,1)) ...
                0.25*(1-IP(:,1)).*(2*IP(:,2)-IP(:,1)) ...
                -0.5*(1-IP(:,1).^2) ...
                -IP(:,2).*(1+IP(:,1)) ...
                0.5*(1-IP(:,1).^2) ...
                -IP(:,2).*(1-IP(:,1))];
        end
        
    otherwise
        error('Only linear and quadratic shape functions are implemented');
end


% Global coordinates of the input points.
xq=psi*elknxyzb(:,1);
yq=psi*elknxyzb(:,2);
zq=psi*elknxyzb(:,3);

if nargout>4
    % Elements of the Jacobian matrix for all input points: dr/de1 and dr/de2, r=(x,y,z).
    dxde1=dNde1*elknxyzb(:,1);
    dyde1=dNde1*elknxyzb(:,2);
    dzde1=dNde1*elknxyzb(:,3);
    
    dxde2=dNde2*elknxyzb(:,1);
    dyde2=dNde2*elknxyzb(:,2);
    dzde2=dNde2*elknxyzb(:,3);
    
    % Normal vector at input points: scalar product (dr/de1 x dr/de2).
    IntExt=sign(elknxyzb(1,4)); % Change sign if the domain is interior. (VC)
    nx=(dyde1.*dzde2-dzde1.*dyde2)*IntExt;
    ny=(dzde1.*dxde2-dxde1.*dzde2)*IntExt;
    nz=(dxde1.*dyde2-dyde1.*dxde2)*IntExt;
end
