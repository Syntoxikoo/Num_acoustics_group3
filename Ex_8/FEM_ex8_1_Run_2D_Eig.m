%==================================================================
% Exercise 9.
%
% Solve a 2D Acoustic Eigenvalue Problem using Linear elements
% 
%==================================================================
clc; close all

% Load Element topology
load('FEM_ex8_1_MeshInfo.mat')
% Plot element
figure(1)
FEM_ex8_1_DrawElementTopology(ex,ey,edof(:,1))

%% Acoustic parameter
c = 343; 
rho = 1.25;
ep = [c rho 2]; %ep =[speed of sound, density, integration rule]

ndof = max(max(edof));
K = zeros(ndof,ndof);
M = zeros(ndof,ndof);

%% Assemble
for i=1:length(ex)
    [ke,me]=FEM_ex8_1_AcoQ4(ex(i,:),ey(i,:),ep); % Generate Element matrix
    % TODO: Use `edof` to add ke and me to the correct places in K and M
    L = zeros(4,ndof);
    for jj = 1:4
        L(jj,edof(i,jj+1)) = 1;
    end
    K = K + L'*ke*L;
    M = M + L'*me * L;

end

%% Solve
[X1,D]=eig(K,M);
for j=1:ndof
    mnorm=sqrt(X1(:,j)'*M*X1(:,j));
    X1(:,j)=X1(:,j)/mnorm;
end

[omega,id]=sort(diag(D)); % Eigenvalue
freq=sqrt(omega)/2/pi;
U=X1(:,id); %Eigenvector
      

%% Plot Modeshape
clc; modnr=1;
figure(1)
for j=1:3 
    for i=1:3
        modnr=modnr+1;
        Ed=FEM_ex8_1_extract(edof,U(:,modnr)); 
        Edabs=abs(Ed); 
        const=max(max(Edabs)); Ed=Ed/const;
        h=fill(ex'+10*1.1*(i-1),ey'-4*1.6*(j-1),Ed');
        MOD= num2str(freq(modnr));
        text(10*1.1*(i-1)+4,-4*1.6*(j-1)-1, MOD)
        set(h,'edgecolor','none')
        hold on, 
    end 
end
axis equal, axis off
