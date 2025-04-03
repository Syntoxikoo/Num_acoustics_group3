function [nvect,tvect]=normals(rzb,topology,see)

% [nvect,tvect]=normals(rzb,topology,see)
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
%         (tvect) the tangentials, same rows as rzb, 
%                 columns rho and z coordinates
%
% If see='y', a figure is plotted with the normals as arrows (quiver).

% Vicente Cutanda Henriquez 05-2007
% Tangentials included, VCH 01-2011


M = size(rzb,1);
[N, ncols]=size(topology);
nknel=ncols-1;
x=linspace(-1,1,nknel);

nvect=[];
for iel=1:N % it is assumed that all bodies have more than one element
   elknrzb(1:nknel,:)=rzb(topology(iel,1:nknel),:);
   [psi, rhoq, zq, nrho, nz]=elemshape(elknrzb,x);
   bbi=find(rzb(:,end)==topology(iel,end));  % find nodes in current body

   if iel==1 | topology(iel,end)~=topology(iel-1,end) % first element in a body
      nvect=[nvect ; nrho(1:end-1,1) nz(1:end-1,1)];
      nold=[nrho(end) nz(end)]/sqrt(nrho(end).^2+nz(end).^2);
      
   elseif iel==N | topology(iel,end)~=topology(iel+1,end) % last element of a body
      nini=nold+[nrho(1) nz(1)]/sqrt(nrho(1).^2+nz(1).^2);
      if ~any(nini)
         disp(['Normal vector not defined at node' num2str(iel)]);
         nini=[NaN NaN];
      end
      if abs(rzb(bbi(1),1))<100*eps & abs(rzb(bbi(end),1))<100*eps % the body includes the z-axis
          nvect=[nvect ; nini ; nrho(2:end,1) nz(2:end,1)];
      elseif abs(rzb(bbi(1),1))>eps & abs(rzb(bbi(end),1))>eps % the body does not include the z-axis (with central hole)
          nvect=[nvect ; nini ; nrho(2:end-1,1) nz(2:end-1,1)];
          nvect(bbi(1),:)=nvect(bbi(1),:)/sqrt(sum(nvect(bbi(1),:).^2))+[nrho(end) nz(end)]/sqrt(nrho(end).^2+nz(end).^2);
      else
          error('Mesh error')
      end
      
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

for kk=1:M % create tangential vector from normal vector
    if rzb(kk,end)>0 % it will point in the clockwise direction
        ph=-pi/2;
    else
        ph=pi/2;
    end
    tvect(kk,1)=nvect(kk,1)*cos(ph)-nvect(kk,2)*sin(ph);
    tvect(kk,2)=nvect(kk,1)*sin(ph)+nvect(kk,2)*cos(ph);
end % Bronshtein, p 230(ES), axis rotation

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
    quiver(rzb(:,1),rzb(:,2),nvect(:,1)*normplot,nvect(:,2)*normplot,0); % plot normal vectors
    quiver(rzb(:,1),rzb(:,2),tvect(:,1)*normplot,tvect(:,2)*normplot,0,'g'); % plot tangential vectors
    hold off;
    grid;
    axis equal;
    title(['Nodes = ' num2str(size(rzb,1)) '  Elements = ' num2str(size(topology,1)) ...
        '  Bodies = ' num2str(num_bodies)]);
    xlabel('rho');
    ylabel('z');
    axis equal
end

