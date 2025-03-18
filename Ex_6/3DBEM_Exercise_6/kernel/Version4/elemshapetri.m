function [psi, xq, yq, zq, nx, ny, nz]=elemshapetri(elknxyzb,IP)

% [psi, xq, yq, zq, nx, ny, nz]=elemshapetri_new(elknxyzb,IP)
%
% Calculate coordinates and normalvector at 'x'
% in an isoparametric triangular element. (linear).
%
%  Input:
%    IP:      real matrix containing the ordinates for the
%             integration points (IP's) x values in 1st column y values in 2nd
%    elknxyzb: real matrix, one row for each node in the element
%             each row contains (x,y,z,body) for the node
%
%  Output:
%    psi:     Real matrix containing one shape function in each
%             column (no. of columns=no. of nodes in an element)
%    xq:      Real column vector containing the x-coordinates
%             of the element for each IP
%    yq:      Real column vector containing the y-coordinates
%             of the element for each IP
%    zq:      Real column vector containing the z-coordinates of
%             the element for each IP
%    nx,ny,nz: 
%             The column vectors containing respectively the x,y
%             and z-coordinates of the normal vector.
%             Note that the normal vector does not have unit length.
%             This is done to avoid division by the jacobian which
%             occasionally is zero. Subsequent multiplication must
%             therefore be omitted.

% Light version with less calculation enabled, to be used in element
% subdivision for near-singular integrals. Vicente Cutanda Henriquez 12-2010.

[nknel,dummy]=size(elknxyzb);
if nknel==3
    % Linear Shape functions:
    psi=[IP(:,1) ...
        IP(:,2) ...
        1-IP(:,1)-IP(:,2)];
    if nargout>4
        % Linear Shape function derivatives
        % Well, they are actually constant
        ettaller=ones(size(IP,1),1);
        nuller=zeros(size(IP,1),1);
        
        dNx=[ettaller nuller -ettaller];
        dNy=[nuller ettaller -ettaller];
        
        %    dNx(:,1)=ones(size(IP,1),1);
        %    dNx(:,2)=zeros(size(IP,1),1);
        %    dNx(:,3)=-ones(size(IP,1),1);
        %
        %    dNy(:,1)=zeros(size(IP,1),1);
        %    dNy(:,2)=ones(size(IP,1),1);
        %    dNy(:,3)=-ones(size(IP,1),1);
    end
elseif nknel==6
    % Quadratic shape functions
    xi1=IP(:,1); xi2=IP(:,2); xi3=1-xi1-xi2;
    nuller=zeros(size(IP,1),1);
    psi=[xi1.*(2*xi1-1) ...
        xi2.*(2*xi2-1) ...
        xi3.*(2*xi3-1) ...
        4*xi1.*xi2 ...
        4*xi2.*xi3 ...
        4*xi1.*xi3];
    if nargout>4
        dNx=[4*xi1-1 ...
            nuller ...
            -(4*xi3-1) ...
            4*xi2 ...
            -4*xi2 ...
            4*(xi3-xi1)];
        dNy=[nuller ...
            4*xi2-1 ...
            -(4*xi3-1) ...
            4*xi1 ...
            4*(xi3-xi2) ...
            -4*xi1];
    end
else
    error('Only linear and quadratic triangular shape functions are implemented');
end

xq=psi*elknxyzb(:,1);
yq=psi*elknxyzb(:,2);
zq=psi*elknxyzb(:,3);

if nargout>4
    % Elements of the Jacobian matrix for all input points: dr/de1 and dr/de2, r=(x,y,z).
    dxde1=dNx*elknxyzb(:,1);
    dyde1=dNx*elknxyzb(:,2);
    dzde1=dNx*elknxyzb(:,3);
    
    dxde2=dNy*elknxyzb(:,1);
    dyde2=dNy*elknxyzb(:,2);
    dzde2=dNy*elknxyzb(:,3);
    
    % Normal vector at input points: scalar product (dr/de1 x dr/de2).
    IntExt=sign(elknxyzb(1,4)); % Change sign if the domain is interior. (VC)
    nx=(dyde1.*dzde2-dzde1.*dyde2)*IntExt;
    ny=(dzde1.*dxde2-dxde1.*dzde2)*IntExt;
    nz=(dxde1.*dyde2-dyde1.*dxde2)*IntExt;
end


% IntExt=sign(elknxyzb(1,4)); % Change direction of the normal vector if the domain is interior. 
% nx=(dNx*elknxyzb(:,2)).*(dNy*elknxyzb(:,3))-(dNy*elknxyzb(:,2)).*(dNx*elknxyzb(:,3))*IntExt;
% ny=(dNy*elknxyzb(:,1)).*(dNx*elknxyzb(:,3))-(dNx*elknxyzb(:,1)).*(dNy*elknxyzb(:,3))*IntExt;
% nz=(dNx*elknxyzb(:,1)).*(dNy*elknxyzb(:,2))-(dNy*elknxyzb(:,1)).*(dNx*elknxyzb(:,2))*IntExt;
