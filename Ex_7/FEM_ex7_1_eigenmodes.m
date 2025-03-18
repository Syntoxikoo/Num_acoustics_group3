%==================================================================
% Exercise 7.1
%
% Solve a 1D Acoustic Eigenvalue Problem using Linear elements
% Boundary conditions: Fixed pressure at both ends%
%
%==================================================================
close all; clc

% Problem data
Lx = 1;     % Length of the domain; fixed arbitrary to 1 (m)
c0 = 342;   % Speed of sound (m/s)

% Step 1: Mesh
ne  = input(' Number of linear elements: ');
nnt = ne+1;     % Total number of nodes
h = Lx/ne;      % Length of the elements
x = 0:h:Lx;     % Coordinates table
nmodes = 6;     % Number of eigensolutions

% Step 2: Compute Elementary matrices
Ke = c0^2* [1,-1;-1,1]*1/h;
Me = [2,1;1,2]*h/6;

% Step 3: Assembling
I = eye(2,2);
K = zeros(nnt,nnt); 
M = zeros(nnt,nnt);
for ie=1: ne
    L = zeros(2,nnt); L(:,ie:ie+1)=I; % Location matrix for element ie
    K = K + L'*Ke*L;
    M = M + L'*Me*L;
end

% Step 4: Boundary conditions: p(1) = p(nnt) = 0
K = K(2:nnt-1,2:nnt-1);
M = M(2:nnt-1,2:nnt-1);

% Step 5: Compute the eigenvalues and eigenvectors
[V,D]   = eig(K,M); 
D       = sqrt(D);
% Nor     = V'*M*V; 
% V       = V*sqrt(inv(Nor)); % Normalization of the eigenvectors
w       = diag(D);
[w,ind] = sort(w); % Sort the eigenvalues in ascending order


% Step 6: Comparison with the exact solution for a selected number of modes
% comparison of the mode shapes
xt = [0:0.025:1] * Lx;
for m=1:nmodes
    figure;
    yt= sqrt(2)*sin(m*pi*xt);
    w_ex = Lx*c0*pi*m;
    plot(xt,yt,'LineWidth',2,'Color',[0 0 0])
    hold on, % Theoretical solution for

    y=[0,real(V(:,ind(m))'),0]; 
    if (yt(2) >0 && y(2) < 0)
        y = -y; 
    end % Recall : p(0)=p(L)=0

    plot(x,y,'-*','LineWidth',2,'Color',[0.5 0.5 0.5]) % FE solution
    title(['Mode' num2str(m)]);
        xlabel('Position x/L');
    ylabel(' Normalized Amplitude ');
    legend('Analytical', 'FEM Linear');
end

%% Conv rate of the algo

Lx = 1;
c0 = 343;
Ne_arr = 12:200;
nmodes = 10;
w_ex = Lx*c0*pi*(1:nmodes);


Deltap = zeros(nmodes, length(Ne_arr));
DeltaWi = zeros(nmodes, length(Ne_arr));
for ll = 1:length(Ne_arr)
    disp("iter :"+ ll)

    ne = Ne_arr(ll);
    nnt=  ne+1;
    h = Lx ./ne;


    p_ex = zeros(nmodes, nnt);
    xt = linspace(0,1,nnt);
    for ii = 1: nmodes
        p_ex(ii,:) = sqrt(2)*sin(ii.*pi.*xt);
    end

    % Step 2: Compute Elementary matrices
    Ke = c0^2* [1,-1;-1,1]*1/h;
    Me = [2,1;1,2]*h/6;
    
    % Step 3: Assembling
    I = eye(2,2);
    K = zeros(nnt,nnt); 
    M = zeros(nnt,nnt);

    for ie=1: ne
        L = zeros(2,nnt); L(:,ie:ie+1)=I; % Location matrix for element ie
        K = K + L'*Ke*L;
        M = M + L'*Me*L;
    end
    % Boundary condition
    K = K(2:nnt-1,2:nnt-1);
    M = M(2:nnt-1,2:nnt-1);

    % Compute eig value & vector
    [V,D]   = eig(K,M); 
    D       = sqrt(D); 
    w       = diag(D);
    [w,ind] = sort(w);
    
    DeltaWi(:,ll) = abs(w(1:10) -w_ex.').^2 ./ abs(w_ex.').^2;

    matched_indices = dsearchn(w_ex', w(1:nmodes));
    for jj =1:nmodes
        y=[0,real(V(:,ind(jj))'),0]; 
        true_mode = matched_indices(jj);
        dot_prod = p_ex(true_mode,2:end-1) * V(:,ind(jj));
        if dot_prod < 0
            y = -y;
        end
    Deltap(jj,ll) = trapz(xt,abs(p_ex(jj,:) - y).^2)./trapz(xt,abs(p_ex(jj,:)).^2);
    end
end

%% Plot convergence

figure; 
hold on;grid on
for m = 1:6
    plot(Ne_arr, DeltaWi(m,:), "LineWidth",1, DisplayName= "mode n°"+m)
end
xlabel("Number of elements")
legend()
set(gca, "Yscale","log")
ylabel("|error|")
xlim([min(Ne_arr) max(Ne_arr)])
title("Convergence of the eigenfrequency in function of the number of element")

figure;
hold on ; grid on 
for m = 1:6
    plot(Ne_arr, Deltap(m,:), "LineWidth",1,DisplayName= "mode n°"+m)
end
xlabel("Number of elements")
ylabel("|error|")
xlim([min(Ne_arr) max(Ne_arr)])
legend()
set(gca, "Yscale","log")
title("Convergence of the eigenfrequency in function of the number of element")


%% fluid particle velocity from Euler eq

c0 = 343;
Lx = 1;
% ne = input("number of element");
ne = 50;
nnt = ne+1;
h = Lx/ne;
x = 0:h:Lx;
nmodes = 6;

% 1D case
Ke = c0^2/h * [1,-1;-1,1] ;
Me = h/6 * [2,1;1,2];

% Define localisation mtx and Assemble
I = eye(2,2);
K = zeros(nnt);
M = zeros(nnt);

for ie = 1:ne
    L = zeros(2, nnt); L(:,ie:ie+1)= I;
    K = K + L'*Ke*L;
    M = M + L'*Me *L;
end

% Boundary
K = K(2:nnt-1,2:nnt-1);
M = M(2:nnt-1,2:nnt-1);

%finding frequency of the modes and pressure 
% response function for each of them
[V,W] = eig(K,M);
W = sqrt(diag(W));
[w, ind] = sort(W);

    

w_ex = (c0*pi/Lx)*(1:nmodes)';
p_ex = sqrt(2)*sin((1:nmodes)'*pi*x/Lx)';

p_n = zeros(nnt, nmodes);
for nn = 1 : nmodes
    p_n(2:end-1,nn) = V(:,ind(nn));
    if dot(p_n(:,nn),p_ex(:,nn).') <0
        p_n(:,nn) = - p_n(:,nn);
    end
end

% velocity 
rho = 1.21;
v_ex = -1j * diff(p_ex)./ (rho .* w_ex.');
v_n = -1j * diff(p_n)./ (rho .* w(1:nmodes)');
xMp = (x(1:end-1) + x(2:end))/2;

m = 6;
figure;
nexttile
plot(xMp,abs(v_n(:,m)), "LineWidth",1, "Color","r"); hold on; grid on
plot(xMp,abs(v_ex(:,m)), "LineWidth",1,"LineStyle","--", "Color","r")

hold off
nexttile
plot(x,abs(p_n(:,m)), "LineWidth",1, "Color","b"); hold on 
plot(x,abs(p_ex(:,m)), "LineWidth",1,"LineStyle","--", "Color","b")
hold off
legend("Vel num","Vel th", "pressure num", "pressure th")


%% Quadratic element
c0 = 343;
Lx = 1;
% ne = input("number of element");


quadK = [7 -8 1; -8 16 -8; 1 -8 7];
quadM = [4 2 -1; 2 16 2; -1 2 4];

%convergence Study
Ne_arr = 12:200;
nmodes = 10;
w_ex = Lx*c0*pi*(1:nmodes);

Deltap = zeros(nmodes, length(Ne_arr));
DeltaWi = zeros(nmodes, length(Ne_arr));

for ll = 1: length(Ne_arr)
    disp("iter :"+ Ne_arr(ll))

    ne = Ne_arr(ll);
    nnt=  ne*2 +1;
    h = Lx ./ne;


    p_ex = zeros(nmodes, nnt);
    xt = linspace(0,1,nnt);
    for ii = 1: nmodes
        p_ex(ii,:) = sqrt(2)*sin(ii.*pi.*xt);
    end

    % Step 2: Compute Elementary matrices
    Ke = c0^2* quadK*1/(3*h);
    Me = quadM*h/30;
    
    % Step 3: Assembling
    K = zeros(nnt,nnt); 
    M = zeros(nnt,nnt);
    for ie = 1:ne
        idx = 2*ie-1;
        
        L = zeros(3, nnt);
        L(:, idx:idx+2) = eye(3);
        
        K = K + L'*Ke*L;
        M = M + L'*Me*L;
    end
    % Boundary condition
    K = K(2:nnt-1,2:nnt-1);
    M = M(2:nnt-1,2:nnt-1);

    % Compute eig value & vector
    [V,D]   = eig(K,M); 
    D       = sqrt(D); 
    w       = diag(D);
    [w,ind] = sort(w);
    
    DeltaWi(:,ll) = abs(w(1:10) -w_ex.').^2 ./ abs(w_ex.').^2;

    matched_indices = dsearchn(w_ex', w(1:nmodes));
    for jj =1:nmodes
        y=[0,real(V(:,ind(jj))'),0]; 
        true_mode = matched_indices(jj);
        dot_prod = p_ex(true_mode,2:end-1) * V(:,ind(jj));
        if dot_prod < 0
            y = -y;
        end
    Deltap(jj,ll) = trapz(xt,abs(p_ex(jj,:) - y).^2)./trapz(xt,abs(p_ex(jj,:)).^2);
    end
end


figure; 
hold on;grid on
for m = 1:6
    plot(Ne_arr, DeltaWi(m,:), "LineWidth",1, DisplayName= "mode n°"+m)
end
xlabel("Number of elements")
legend()
set(gca, "Yscale","log")
ylabel("|error|")
xlim([min(Ne_arr) max(Ne_arr)])
title("Convergence of the eigenfrequency in function of the number of element")

figure;
hold on ; grid on 
for m = 1:6
    plot(Ne_arr, Deltap(m,:), "LineWidth",1,DisplayName= "mode n°"+m)
end
xlabel("Number of elements")
ylabel("|error|")
xlim([min(Ne_arr) max(Ne_arr)])
legend()
set(gca, "Yscale","log")
title("Convergence of the eigenfrequency in function of the number of element")
hold off



ne = 50;
nnt = ne*2+1;
h = Lx/ne;
x = 0:h/2:Lx;
nmodes = 6;
% 1D caseL
Ke = c0^2/h * quadK;
Me = h/6 * quadM;

% Define localisation mtx and Assemble
K = zeros(nnt);
M = zeros(nnt);

for ie = 1:ne
    idx = 2*ie-1;
    
    L = zeros(3, nnt);
    L(:, idx:idx+2) = eye(3);
    
    K = K + L'*Ke*L;
    M = M + L'*Me*L;
end

% Boundary
K = K(2:nnt-1,2:nnt-1);
M = M(2:nnt-1,2:nnt-1);

%finding frequency of the modes and pressure 
% response function for each of them
[V,W] = eig(K,M);
Nor     = V'*M*V; 
V       = V*sqrt(inv(Nor));
W = sqrt(diag(W));
[w, ind] = sort(W);

    

w_ex = (c0*pi/Lx)*(1:nmodes)';
p_ex = sqrt(2)*sin((1:nmodes)'*pi*x/Lx)';

p_n = zeros(nnt, nmodes);
for nn = 1 : nmodes
    p_n(2:end-1,nn) = V(:,ind(nn));
    if dot(p_n(:,nn),p_ex(:,nn).') <0
        p_n(:,nn) = - p_n(:,nn);
    end
end

% velocity 
rho = 1.21;
v_ex = -1j * diff(p_ex)./ (rho .* w_ex.');
v_n = -1j * diff(p_n)./ (rho .* w(1:nmodes)');
xMp = (x(1:end-1) + x(2:end))/2;

m = 1;
figure;
nexttile
plot(xMp,abs(v_n(:,m)), "LineWidth",1, "Color","r"); hold on; grid on
plot(xMp,abs(v_ex(:,m)), "LineWidth",1,"LineStyle","--", "Color","r")

hold off
nexttile
plot(x,abs(p_n(:,m)), "LineWidth",1, "Color","b"); hold on 
plot(x,abs(p_ex(:,m)), "LineWidth",1,"LineStyle","--", "Color","b")
hold off
legend("Vel num","Vel th", "pressure num", "pressure th")

