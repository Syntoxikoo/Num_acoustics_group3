function []=plotgeom(rzb,topology)

% []=plotgeom(rzb,topology)
%
% Creates a figure of the generator of an axisymmetrical mesh
%
% Input variables:
%     -rzb:      node positions, first column is the rho-coordinate, second
%                column is z-coordinate, and third column is the body
%                number to which the node belongs to.
%     -topology: each row contains the node numbers (row number in rzb) of the
%                nodes in one element. The last column is the body number the
%                element belongs to.

% Vicente Cutanda Henríquez, 07-2011
 

% Plot geometry
%    figure;
    for bb=1:max(abs(topology(:,end)))
        eleb=find(abs(topology(:,end))==bb);
        nodb=find(abs(rzb(:,end))==bb);
        plot(rzb(nodb,1),rzb(nodb,2),'r:');
        hold on;
        for nn=1:length(nodb)
            if any(topology(eleb,1)==nodb(nn))|any(topology(eleb,3)==nodb(nn))
                plot(rzb(nodb(nn),1),rzb(nodb(nn),2),'ko');
            else
                plot(rzb(nodb(nn),1),rzb(nodb(nn),2),'k+');
            end
        end
    end
    hold off;
    grid;
    axis equal;
%     title(['Nodes = ' num2str(size(rzb,1)) '  Elements = ' num2str(size(topology,1)) ...
%         '  Bodies = ' num2str(num_bodies)]);
    xlabel('rho');
    ylabel('z');
end