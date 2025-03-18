% Calculates the pressure on a vibrating torus:
clear

R1=0.5;              % Inner radius of the torus
R2=1;              % Outer radius of the torus
u0=1;              % Maximum velocity amplitude
k=1;               % Wavenumber, m-1
m=0;               % Excitación ejesimétrica, modo circunferencial m=0
quadelem=1;        % If 0, linear elements (2 nodes) are used instead of quadratic (3 nodes)

% CONDICIONES AMBIENTALES
pa = 101325;          % Presion Atmosferica (Pa)
t = 20;               % Temperatura (ºC)
Hr = 50;              % Humedad relativa (%)
[rho,c,cf,CpCv,nu,alfa]=amb2prop(pa,t,Hr,1000); 

 

% GENERACIÓN DE LA GEOMETRIA DEL DOMINIO
segments=[R1 0 R2 0 10 (R2-R1)/2+eps 0;R2 0 R1 0 10 (R2-R1)/2+eps 0;];
[rzb,topology]=nodegen(segments,'n',{},quadelem); % nodes and elements

% topology(end,3)=topology(1,1);
% rzb=rzb(1:end-1,:);
M=size(rzb,1);N=size(topology,1);     % M nodos, N elementos

figure; plot(rzb(:,1),rzb(:,2),'ro');grid
% Velocity (first-order oscillating sphere)
vf=u0*rzb(:,2);



% <<<<<<<<<<<<< Correr desde el inicio hasta aquí para ver la geometría sin calcular


% CÁLCULO DE LAS MATRICES BEM
% Calcula matrices BEM
[A,B,CConst]=BEMEquat0(rzb,topology,k,m);
disp(['Números de condición, de A: ' num2str(cond(A)) ' y de B: ' num2str(cond(B))])
B=i*k*rho*c*B;
% CALCULO DE LA SOLUCIÓN SOBRE LA SUPERFICIE
pf=A\(-B*vf); % the two test cases


% Plot solution (first-order oscillating sphere)
figure;
plot(1:M,abs(pf),'-rx'); grid
xlabel('Angle');ylabel('|pressure on the surface|');



