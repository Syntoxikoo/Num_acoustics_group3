function [rzb,topology,rzline,segrzb]=nodegen(segments,see,varargin)

% Meshing function in 2D and axisymmetrical BEM.
%
% Basic usage:
% ============
% [xyb,topology]=nodegen(segments,see);
% [rzb,topology]=nodegen(segments,see);
%
% In its simplest form, nodegen builds geometry arrays from a description
% of the geometry as circular or straight segments of constant mesh
% density. The mesh generated consists of quadratic three-node elements. If
% there are different bodies, they are detected and numbered.
%
% If the formulation is axisymmetrical, the mesh represents the generator
% and cylindrical coordinates (rho,z) are assumed. The objects can be
% connected to the z-axis or detached from it as closed toroidal
% geometries.
%
% If the formulation is 2-dimensional, the mesh represents the 2D boundary
% and 2D cartesian coordinates (x,y) are assumed. Objects must then be
% closed.
%
% Input variables:
%     -segments: array with a segment description per row. The segments
%                must be consecutive in each body and in a clockwise
%                direction. The columns are: rho and z of the starting
%                point, rho and z of the end point (x and y in 2D), number
%                of elements in the segment, curvature radius and mesh
%                density (elements/m). 
%                The sign of the radius indicates whether it bends outwards
%                (+) or inwards (-). If the radius is Inf (or 0), a
%                straight segment is used.
%                The necessary number of elements to achieve the density
%                given in column 7 is computed. If this number is larger
%                than the fixed number in column 5, it is used instead.
%     -see:      if is set to 'y' or 'Y', a plot is made with the geometry.
%
% Output variables:
%     -rzb:      node positions, first column is the rho-coordinate, second
%                column is z-coordinate, and third column is the body
%                number to which the node belongs to.
%     -topology: each row contains the node numbers (row number in rzb) of
%                the nodes in one element. The last column is the body number
%                the element belongs to.
% 
% Advanced usage:
% ==================
% [rzb,topology,rzline,segrzb]=nodegen(segments,see,segmfun,quad);
%
% Other features have also been implemented. In addition to the variables
% above, there are new input and output variables:
%
% Additional input variables:
%     -segmfun:  cell vector with the inline functions for segment shapes.
%                Must be a cell, even for a single inline or anonimous
%                function. If the radius (6th column in "segments") is
%                imaginary as (+/-N*j), with N=1,2,3..., then the segment
%                follows the inline (or anonimous) function segmfun{N} . If
%                the sign is +, it is z=z(r), and if the sign is -, it is
%                r=r(z). Care must be taken for the initial and final
%                points to fulfill the function, and for the functions to
%                be single valued within the interval.
%     -quad:     If 1, quadratic elements are used (3 nodes, default). If 0,
%                linear elements (2 nodes) are produced.
%
% Additional output variables:
%     -rzline:   matrix with two columns: 
%                -positions of the nodes in a one-dimensional axis. Equivalent
%                 to extending the generator over a straight line.
%                -corresponding body number
%                If the object is closed, the rzline value of its first 
%                node is not zero, but the distance to its last node.
%     -segrzb:   vector with one number per node, as in rzb, indicating the
%                segment number "m" each node belongs to. Useful for
%                assignment of boundary conditions. If a node belongs to
%                two segments with numbers m and m+1, it is given the
%                number m+1/2 (except for the first node in a closed body).
%
% Additional features:
%     -Varying element length: This feature is triggered by making the
%                number of elements in column 5 of "segments" imaginary.
%                The elements will then vary in length along the segment
%                following a geometrical series, and will match the element
%                sizes of the neighboring segments, or else the default
%                mesh density if there is no neighboring segment. This
%                feature is only implemented for straight segments.
%     -Bézier segment: the four coordinates in the first four columns of 
%                "segments" become the four points of a cubic Bézier curve,
%                p0, p1, p2 and p3 if they are given as complex numbers with
%                real and imaginary parts corresponding to rho and z
%                components. References:
%                *Numerical Recipes:Prautzsch, H., Boehm, W., and Paluszny, M. 2002
%                *Bézier and B-Spline Techniques (Berlin:Springer)

% -Vicente Cutanda Henríquez 1995 
% -Translation from Pascal to Matlab: Johan Gramtorp 2000
% -Vicente Cutanda Henríquez: help text 2004, inline function segments
% 2007, rzline 2010.
% Vicente Cutanda Henríquez, rzline addition, 04-2011
% Vicente Cutanda Henríquez, varying element length, Bézier segments, 06-2011
% Vicente Cutanda Henríquez, objects not attached to the z-axis, 05-2014
% Vicente Cutanda Henríquez, Output corresponding segment numbers. Common
% function in 2D and Axi BEM. 10-2014.
 

if nargin>2
    segmfun=varargin{1};
end

number_seg=size(segments,1);	% number of segments

segments(isinf(segments(:,6)),6)=0; % Allow Inf as an indicator of straight segment

segmentscell=cell(number_seg,3); % stores partial data for every segment: rho, z and length

% Find straigth segments that will have varying element length
m_re=find(imag(segments(:,5))==0 | segments(:,6)~=0);
m_im=find(imag(segments(:,5))~=0 & segments(:,6)==0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=m_re';	% For every segment with uniform element length
    rho_start=segments(m,1);
    z_start=segments(m,2);
    rho_end=segments(m,3);
    z_end=segments(m,4);
    number_elemMIN=segments(m,5);	% number of elements for the actual segment
    radius=segments(m,6);			% radius=0 for straight segments
    meshdens=segments(m,7);			% element density in elements per meter
    
    if ~isreal([rho_start z_start rho_end z_end]) % Bézier curve with extra points as imaginary parts
        % Bézier curves:
        % Ref from Numerical Recipes:Prautzsch, H., Boehm, W., and Paluszny, M. 2002
        % Bézier and B-Spline Techniques (Berlin:Springer)

        p0=rho_start;p1=z_start;p2=rho_end;p3=z_end; % in this case the four coordinates become four points
        t=linspace(0,1,1e6+1)'; % number of function values to evaluate segment length

        % Cubic Bézier curve, parametric form
        B=(1-t).^3*p0+3*(1-t).^2.*t*p1+3*(1-t).*t.^2*p2+t.^3*p3;

        table=[real(B) imag(B)]; % rho and z coordinates
        table=[table [0;sqrt(sum([diff(table).^2]'))']]; % calculate distances
        table(:,4)=cumsum(table(:,3));
        
        % calculate segment length and number of elements
        segleng=sum(table(:,3)); % approximate length of the segment
        number_elem=max(abs(number_elemMIN),ceil(segleng*meshdens));
        
%         step1=segleng/(number_elem); % step, size of the first element
%         [dum,it]=min(abs(table(:,4)-step1));step=t(it); % step in the parameter
%         tt=linspace(0,1,round(1/step));
        tt=linspace(0,1,number_elem);
        
        for pp=1:length(tt)-1
            [dum,ai1]=min(abs(t-tt(pp)));[dum,ai2]=min(abs(t-tt(pp+1)));
            segmentscell{m,1}(pp*2-1,1)=table(ai1,1); % rho 
            segmentscell{m,2}(pp*2-1,1)=table(ai1,2); % z
            segmentscell{m,3}(pp*2-1,1)=table(ai1,4); % length
            % mid nodes must be in the middle:
            [dum,ai3]=min(abs(table(:,4)-(table(ai1,4)+table(ai2,4))/2));
            segmentscell{m,1}(pp*2,1)=table(ai3,1); % rho 
            segmentscell{m,2}(pp*2,1)=table(ai3,2); % z
            segmentscell{m,3}(pp*2,1)=table(ai3,4); % length
        end
        % last node:
        segmentscell{m,1}(length(tt)*2-1,1)=table(end,1); % rho
        segmentscell{m,2}(length(tt)*2-1,1)=table(end,2); % z
        segmentscell{m,3}(length(tt)*2-1,1)=table(end,4); % length
        
        
    elseif radius==0; % straigth segment
        % calculate segment length and number of elements
        segleng=sqrt((rho_end-rho_start)^2+(z_end-z_start)^2);
        number_elem=max(abs(number_elemMIN),ceil(segleng*meshdens));
        
        % calculate the rho's
        segmentscell{m,1}=linspace(rho_start,rho_end,number_elem*2+1)';
        % calculate the z's
        segmentscell{m,2}=linspace(z_start,z_end,number_elem*2+1)';
        
        % Calculate element lengths
        step=sqrt((rho_end-rho_start)^2+(z_end-z_start)^2)/(number_elem*2);
        segmentscell{m,3}=(0:number_elem*2)'*step;
        
    elseif isreal(radius) % Radius <> 0 and real
        distance=sqrt((rho_end-rho_start)^2+(z_end-z_start)^2);
        diameter=abs(2*radius);
        if distance>abs(2*radius);	% circle not big enough !!!
            error('Wrong input data: Circle not big enough')
            distance
            diameter
        else
            % get the coordinates of the 2 circumference centres
            % starting point and end point has the same z-value
            if z_end==z_start;
                % find the two centers
                rho01=(rho_end-rho_start)/2+rho_start;
                rho02=rho01;
                z01=z_start-sqrt(radius^2-((rho_end-rho_start)/2)^2);
                z02=z_start+sqrt(radius^2-((rho_end-rho_start)/2)^2);
            else	% different z-values
                % find the two centers
                m0=-(rho_end-rho_start)/(z_end-z_start);
                b0=(z_start+z_end)/2-m0*(rho_start+rho_end)/2;
                a=1+m0^2;
                b=2*m0*(b0-z_start)-2*rho_start;
                c=rho_start^2+(b0-z_start)^2-radius^2;
                rho01=(-b+sqrt(b^2-4*a*c))/(2*a);
                rho02=(-b-sqrt(b^2-4*a*c))/(2*a);
                z01=m0*rho01+b0;
                z02=m0*rho02+b0;
            end
%             % decide which curvature centre will be used
%             ang=angle(rho_end-rho_start+i*(z_end-z_start));
%             if ang<0;
%                 ang=2*pi+ang;
%             end
%             if ang>pi/2;
%                 cond=0;
%             else
%                 cond=1;
%             end
%             if (cond&(radius>0))|((~cond)&(radius<0));
%                 rho0=rho01;
%                 z0=z01;
%             else
%                 rho0=rho02;
%                 z0=z02;
%             end
            % decide which curvature centre will be used
            vecseg=[rho_end-rho_start z_end-z_start 0]; % Sometimes it does not work. Can changed back to the previous version above.
            vec1=[rho01-rho_start z01-z_start 0];
            vec2=[rho02-rho_start z02-z_start 0];
            cvec1=cross(vecseg,vec1);cvec1=(sign(cvec1(end))+1)/2;
            cvec2=cross(vecseg,vec2);cvec2=(sign(cvec2(end))+1)/2;
            if cvec1==cvec2, error('Geometrical error arc segment');end 
            
            if (cvec1&&(radius>0))||((~cvec1)&&(radius<0))
                rho0=rho02;
                z0=z02;
            else
                rho0=rho01;
                z0=z01;
            end
            % change to polar coordinates centered in the circumference
            r_start=sqrt((rho_start-rho0)^2+(z_start-z0)^2);
            r_end=sqrt((rho_end-rho0)^2+(z_end-z0)^2);
            theta_start=angle(rho_start-rho0+i*(z_start-z0));
            theta_end=angle(rho_end-rho0+i*(z_end-z0));
            % calculate node positions
            ang=theta_end-theta_start;
            if ang<0;
                ang=2*pi+ang;
            end
            if ang>pi;
                signo=-1;
                ang=2*pi-ang;
            else
                signo=1;
            end
            if (diameter==distance)&&(radius>0)
                signo=-1;
            end
            if (diameter==distance)&&(radius<0)
                signo=1;
            end
            
            % calculate segment length and number of elements
            segleng=ang*r_start;
            number_elem=max(abs(number_elemMIN),ceil(segleng*meshdens));
            
            step=ang/(number_elem*2);
            h=theta_start+signo*(0:number_elem*2)'*step;
%            segmentscell{m,1}=abs(r_start*cos(h)+rho0);
            segmentscell{m,1}=r_start*cos(h)+rho0;
            segmentscell{m,2}=r_start*sin(h)+z0;
            step=abs(radius*ang)/(number_elem*2);
            segmentscell{m,3}=(0:number_elem*2)'*step;
        end
        
    else % remaining case: radius imaginary and non-zero (01/2007 addition)
        
        sfun=segmfun{abs(round(imag(radius)))};
        n_ele=10000; % number of function values to evaluate segment length
        if sign(imag(radius))==1 % table of rho-z values for many points
            rr=linspace(rho_start,rho_end,n_ele);
            table=[rr' sfun(rr)'];
        else
            zz=linspace(z_start,z_end,n_ele);
            table=[sfun(zz)' zz'];
        end
        
        table=[table sqrt(sum([0 0 ; (diff(table)).^2]'))']; % calculate distances
        for pp=1:size(table,1); table(pp,4)=sum(table(1:pp,3));end
        
        % calculate segment length and number of elements
        segleng=sum(table(:,3)); % approximate length of the segment
        number_elem=max(abs(number_elemMIN),ceil(segleng*meshdens));
        
        step=segleng/(number_elem*2); % step distance between nodes
        for pp=0:number_elem*2
            [aa,ai]=min(abs(table(:,4)-pp*step));
            segmentscell{m,1}(pp+1,1)=table(ai,1);
            segmentscell{m,2}(pp+1,1)=table(ai,2);
        end
        segmentscell{m,3}=(0:number_elem*2)'*step;
        
    end
end % m loop, uniform segment length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=m_im';	% For every segment with varying element length
    rho_start=real(segments(m,1));
    z_start=real(segments(m,2));
    rho_end=real(segments(m,3));
    z_end=real(segments(m,4));
    number_elemMIN=segments(m,5);	% number of elements for the actual segment
    meshdens=segments(m,7);			% element density in elements per meter
    
    % Obtain the reference element sizes (la, lz) from the neighbouring segments, if they exist.
    la=0; lz=0;
    for mm=1:length(m_re)
        if abs(rho_start - segmentscell{m_re(mm),1}(end)) < eps && abs(z_start - segmentscell{m_re(mm),2}(end)) < eps
            la=segmentscell{m_re(mm),3}(end)-segmentscell{m_re(mm),3}(end-2);
        end
        if abs(rho_end - segmentscell{m_re(mm),1}(1)) < eps && abs(z_end - segmentscell{m_re(mm),2}(1)) < eps
            lz=segmentscell{m_re(mm),3}(3)-segmentscell{m_re(mm),3}(1);
        end
    end
    
    % Compute straight segments. Other kinds of segments not yet implemented for varying element length:

    % calculate segment length and number of elements
    segleng=sqrt((rho_end-rho_start)^2+(z_end-z_start)^2);
    number_elem=max(abs(number_elemMIN),ceil(segleng*meshdens));
    lm=segleng/number_elem; % length of one element if all elements are the same size in the segment
    
    
    if (la==0 && lz==0) || abs(la-lz)<eps % if there are no reference neighboring segments or the references are equal, division is uniform
        % calculate the rho's
        segmentscell{m,1}=linspace(rho_start,rho_end,number_elem*2+1)';
        % calculate the z's
        segmentscell{m,2}=linspace(z_start,z_end,number_elem*2+1)';
        % calculate element lengths
        step=sqrt((rho_end-rho_start)^2+(z_end-z_start)^2)/(number_elem*2);
        segmentscell{m,3}=(0:number_elem*2)'*step;
    else % varying element length
        if la==0, la=lm;end % if one end has no reference, the default mesh density is used
        if lz==0, lz=lm;end
        
        % create finite geometrical series that matches the segment and neighboring meshes
        ii=round(1/(1-log10(la/lz-(la/lz-1)/(1+segleng/la))/log10(la/lz)));
        pp=zeros(1,ii+1); % order of the polynomial + 1
        pp(1)=segleng;pp(2)=-segleng-la;pp(end)=la;
        rr=roots(pp);
        if la>lz
            kk=rr(imag(rr)==0 & real(rr)>1+1000*eps); kk=kk(1);
        else
            kk=rr(imag(rr)==0 & real(rr)<1-1000*eps); kk=kk(1);
        end
        tmp=cumsum([0 la./kk.^(1:ii-1)])';steps=[];
        for hh=1:length(tmp)-1 % introduce midnodes (quadratic elements assumed)
            steps=[steps tmp(hh) (tmp(hh)+tmp(hh+1))/2];
        end
        steps= [steps tmp(end)]'; % varying positions along the segment
        
        % calculate the rho's
        segmentscell{m,1}=rho_start+steps*((rho_end-rho_start))/segleng;
        % calculate the z's
        segmentscell{m,2}=z_start+steps*((z_end-z_start))/segleng;
        % element lengths
        segmentscell{m,3}=steps;
    end
    
end % m loop, varying element length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Fill the nodes array rzb and the length rzline with the calculated values in the cell segmentscell
body=1;rzb=[];closedbody=0;segrzb=[];
for m=1:number_seg;	% For every segment

    if m==1  % The very first segement of the geometry
        bodyfirst=[segmentscell{1,1}(1) segmentscell{1,2}(1) size(rzb,1)+1];
        rzb=[segmentscell{1,1} segmentscell{1,2} ones(length(segmentscell{1,2}),1)];
        rzline=[segmentscell{1,3} ones(length(segmentscell{1,2}),1)];
        segrzb=0.5; % The very fist node will have 1/2 as segment number. Note that this node also belongs to the last segment in closed bodies.
        segrzb=[segrzb; m*ones(length(segmentscell{1,2})-1,1)]; % Fill vector with the segments every node belongs to
    elseif abs(segmentscell{m,1}(1) - segmentscell{m-1,1}(end)) < 100*eps && abs(segmentscell{m,2}(1) - segmentscell{m-1,2}(end)) < 100*eps % detects if the new segment continues right after the last one
        rzb=[rzb; segmentscell{m,1}(2:end) segmentscell{m,2}(2:end) body*ones(length(segmentscell{m,2})-1,1)];
        rzline=[rzline ; rzline(end,1)+segmentscell{m,3}(2:end) body*ones(length(segmentscell{m,2})-1,1)];
        segrzb(end,1)=segrzb(end,1)+0.5; % The node sharing two segments will have the m-1/2 as segment number (midway between the two numbers)
        segrzb=[segrzb; m*ones(length(segmentscell{m,2})-1,1)]; % Fill vector with the segments every node belongs to
    else % if not, we are starting a new body
        body=body+1; closedbody=[closedbody 0]; 
        bodyfirst=[segmentscell{m,1}(1) segmentscell{m,2}(1) size(rzb,1)+1];
        rzb=[rzb; segmentscell{m,1} segmentscell{m,2} body*ones(length(segmentscell{m,2}),1)];
        rzline=[rzline ; segmentscell{m,3} body*ones(length(segmentscell{m,2}),1)];
        segrzb=[segrzb; 0.5]; % The fist node in the body will have 1/2 as segment number. Note that this node also belongs to the last segment in closed bodies.
        segrzb=[segrzb; m*ones(length(segmentscell{m,2})-1,1)]; % Fill vector with the segments every node belongs to
    end
    
    if abs(segmentscell{m,1}(end) - bodyfirst(1)) < 100*eps && abs(segmentscell{m,2}(end) - bodyfirst(2)) < 100*eps % Bodies away from z-axis
        rzb=rzb(1:end-1,:); % remove last node, which is the same as the first one in the body
        % The first node of a cosed body will carry the segment number of the first segment, not the last !!!
        segrzb=segrzb(1:end-1,:); % remove last node, also from the segments list.
        rzline(bodyfirst(3),1)=rzline(end,1); % First body node will hold the last position, instead of a zero. This will be handled in susequent functions.
        rzline=rzline(1:end-1,:);
        closedbody(body)=1;
    end
end


% Calculate number of elements in each body and fill topology matrix
topology=[];
num_bodies=max(rzb(:,3));
for ii=1:num_bodies
    ind=find(rzb(:,3)==ii);
    if closedbody(ii)==1 % Closed bodies
        topology=[topology ; ind(1:2:end-1) ind(2:2:end) [ind(3:2:end) ; ind(1)] ii*ones((length(ind))/2,1)];
    else % Bodies attached to the z-axis
        topology=[topology ; ind(1:2:end-1) ind(2:2:end) ind(3:2:end) ii*ones((length(ind)-1)/2,1)];
    end
end

% Change from quadratic to linear elements
if nargin>3 && varargin{2}~=1
    tmp=zeros(size(topology,1)*2,3);tmp(1:2:end,:)=topology(:,[1 2 4]);tmp(2:2:end,:)=topology(:,[2 3 4]);
    topology=tmp;
    quad=0;
else
    quad=1;
end

% Plot geometry
if see(1)=='y' || see(1)=='Y'
    figure;
    for tt=1:size(topology,1) % Plot segments alternating 3 colors
        if quad
            switch mod(segrzb(topology(tt,2))-1,5)+1
                case 1, colel='r'; case 2, colel='b'; case 3, colel='g'; case 4, colel='m'; case 5, colel='b';
            end
            plot(rzb(topology(tt,1:3),1),rzb(topology(tt,1:3),2),[colel ':'],"LineWidth",1);
            hold on;
            plot(rzb(topology(tt,[1 3]),1),rzb(topology(tt,[1 3]),2),[colel 'o'],"LineWidth",1);
            plot(rzb(topology(tt,2),1),rzb(topology(tt,2),2),[colel '+'],"LineWidth",1);
        else
            switch mod(round(segrzb(topology(tt,1))+0.2)-1,3)+1 % the first node could be shared with another segment, see line 384
                case 1, colel='r'; case 2, colel='b'; case 3, colel='g'; case 4, colel='y'; case 5, colel='m';
            end
            plot(rzb(topology(tt,1:2),1),rzb(topology(tt,1:2),2),[colel ':o'],"LineWidth",1);
            hold on;
%            plot(rzb(topology(tt,1:2),1),rzb(topology(tt,1:2),2),[colel 'o']'ko');
        end
            
    end
    
    hold off;
    grid;
    axis equal;
    title(['Nodes = ' num2str(size(rzb,1)) '  Elements = ' num2str(size(topology,1)) ...
        '  Bodies = ' num2str(num_bodies) '  Segments = ' num2str(number_seg)]);
    xlabel('rho');
    ylabel('z');
end
