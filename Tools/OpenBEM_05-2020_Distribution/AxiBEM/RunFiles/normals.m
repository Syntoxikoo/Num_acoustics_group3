function [nvect]=normals(rzb,topology,see)

% [nvect]=normals(rzb,topology,see)
%
% Calculate normalvectors at the nodes using shape functions,
% element by element. Nodes belonging to two elements get the
% bisecting vector of their two normals.
% The vectors have unit length.
%
% Input: (rzb, topology) the geometry matrices with body numbers
%
% Output: (nvect) the normals, same rows as rzb, 
%                 columns rho and z coordinates
%
% If see='y', a figure is plotted with the normals as arrows (quiver).

% Vicente Cutanda Henriquez 05-2007


M = size(rzb,1);
[N, ncols]=size(topology);
nknel=ncols-1;
x=linspace(-1,1,nknel);

nvect=[];
for iel=1:N % it is assumed that all bodies have more than one element
   elknrzb(1:nknel,:)=rzb(topology(iel,1:nknel),:);
   [psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
   
   if iel==1 % first element in first body
      nvect=[nvect ; nrho(1:end-1,1) nz(1:end-1,1)];
      nold=[nrho(end) nz(end)]/sqrt(nrho(end).^2+nz(end).^2);
      
   elseif topology(iel,end)~=topology(iel-1,end) % first element in other bodies
      nvect=[nvect ; nrho(1:end-1,1) nz(1:end-1,1)];
      nold=[nrho(end) nz(end)]/sqrt(nrho(end).^2+nz(end).^2);
      
   elseif iel==N % last element of last body
      nini=nold+[nrho(1) nz(1)]/sqrt(nrho(1).^2+nz(1).^2);
      if ~any(nini)
         disp(['Normal vector not defined at node' num2str(iel)]);
         nini=[NaN NaN];
      end
      nvect=[nvect ; nini ; nrho(2:end,1) nz(2:end,1)];
      
   elseif topology(iel,end)~=topology(iel+1,end) % last element of other bodies
      nini=nold+[nrho(1) nz(1)]/sqrt(nrho(1).^2+nz(1).^2);
      if ~any(nini)
         disp(['Normal vector not defined at node' num2str(iel)]);
         nini=[NaN NaN];
      end
      nvect=[nvect ; nini ; nrho(2:end,1) nz(2:end,1)];
      
   else % all other elements
      nini=nold+[nrho(1) nz(1)]/sqrt(nrho(1).^2+nz(1).^2);
      if ~any(nini)
         disp(['Normal vector not defined at node' num2str(iel)]);
         nini=[NaN NaN];
      end
      nvect=[nvect ; nini ; nrho(2:end-1,1) nz(2:end-1,1)];
      nold=[nrho(end) nz(end)]/sqrt(nrho(end).^2+nz(end).^2);
   end
   
end

% normalize the vectors
nvect = nvect./(sqrt(nvect(:,1).^2 + nvect(:,2).^2)*[1 1]);

if see(1)=='y' | see(1)=='Y'
    figure;
    num_bodies=abs(topology(end,end));
    for bb=1:num_bodies
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
    normplot=mean(mean(abs(diff(rzb(:,1:2)))));
    quiver(rzb(:,1),rzb(:,2),nvect(:,1)*normplot,nvect(:,2)*normplot,0);
    hold off;
    grid;
    axis equal;
    title(['Nodes = ' num2str(size(rzb,1)) '  Elements = ' num2str(size(topology,1)) ...
        '  Bodies = ' num2str(num_bodies)]);
    xlabel('rho');
    ylabel('z');
    axis equal
end

