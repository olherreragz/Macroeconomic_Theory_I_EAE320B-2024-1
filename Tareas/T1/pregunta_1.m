%% Pregunta 1 Tarea 1 Teoria Macroeconomica
clc; clear; close all;
tx= {'Interpreter','Latex','FontSize', 15};
% PARAMETROS
A = 1.25;
pI = 1;
F = 0.475;
gamma = 0.075;
beta = 0.95;
delta = 0.08;
alpha = 0.48;
kmin = 0.05;
kmax = 20;
Nk = 1000;
max_iter = 1000;
tol = 1e-6; 
k = linspace(kmin,kmax,Nk);


tic
[Cti_L,kti_L,Vf_inv,kp_L, iter_L] = value_inversion(beta, alpha, delta, A, max_iter, k, tol,pI,F,gamma);
toc

tic
[Vf_no, policy_choice, iter] = value_function_sin_inversion(beta, alpha, delta, A, max_iter, k, tol);
toc 
% Ahora computamos la funcion de valor mixta, quem seria tomar los valores
% mas altos de las 2
Vf_mixta = max(Vf_no, Vf_inv');
figure;
hold on;
box on 
plot(k,Vf_no, "r")
plot(k,Vf_inv, "b")
plot(k, Vf_mixta, 'LineWidth', 8, 'Color', [0.75 1 0.75 0.5]);
title('Value Functions', tx{:})
xlabel('$Capital\ K_t$', tx{:})
ylabel('$Valor\ Value\ Function$', tx{:})
legend('Sin Inversion $V^i(K_t)$','Con Inversion $V^a(K_t)$','Con Opcion de Inversion $V(K_t)$','Location', 'best', tx{:}, 'FontSize',13);
hold off
%% 
%Creamos vector indicador
vector_indicador = ones(size(Vf_no));
vector_indicador(Vf_inv' > Vf_no) = 1;
vector_indicador(~(Vf_inv' > Vf_no)) = 0;


i_inv = k(kp_L) - (1-delta)*k;
i_inv(vector_indicador == 0) = 0 ;

% Grafique la policy function para esta planta
figure;
box on;
plot(k,i_inv,"r")
title('Policy Function de la Inversion', tx{:})
xlabel('Capital $K_t$', tx{:})
ylabel('Inversion $I_t$', tx{:})
% legend('Inversion','Location', 'best', tx{:});

%% Pregunta 2
% 2.1
delta1 = 0.08;
delta2 = 0.25;
pdelta1 = 0.7;
pdelta2 = 0.3;
[Cti,kti,Vf_inv_estocastico,kp, iter] = value_inversion_estocastica(beta, alpha, delta1,delta2, A, max_iter, k, tol,pI,F,gamma, pdelta1,pdelta2);
[Vf_no_inv_estocastico, policy_choice, iter1] = value_function_sin_inversion_estocastica(beta, alpha, delta1,delta2,pdelta1,pdelta2, A, max_iter, k, tol);
figure;
hold on;
box on;
plot(k, Vf_inv_estocastico, "b")
plot(k,Vf_no_inv_estocastico, "r")
Vf_mixta_estocastica = max(Vf_no_inv_estocastico, Vf_inv_estocastico');
plot(k, Vf_mixta_estocastica, 'LineWidth', 8, 'Color', [0.75 1 0.75 0.5]);
title('Value Functions Estocasticas', tx{:})
xlabel('$Capital\ K_t$', tx{:})
ylabel('$Valor\ Value\ Function$', tx{:})
legend(...
    'V.F. Estocastica Con Inversion $E_0 (V^a(K_t))$',...
    'V.F. Estocastica Sin Inversion $E_0 (V^i(K_t))$',...
    'V.F. Estocastica Con Opcion de Inversion $E_0(V(K_t))$',...
    'Location','best',...
    tx{:},...
    'FontSize',10 ...
        );
hold off;

% %Creamos vector indicador para el caso estocástico
% vector_indicador_estocastico = ones(size(Vf_no_inv_estocastico));
% vector_indicador_estocastico(Vf_inv_estocastico' > Vf_no_inv_estocastico) = 1;
% vector_indicador_estocastico(~(Vf_inv_estocastico' > Vf_no_inv_estocastico)) = 0;
% 
% expected_delta = pdelta1*delta1 + pdelta2*delta2;
% 
% i_inv_estocastico = k(kp) - (1-expected_delta)*k;
% i_inv_estocastico(vector_indicador_estocastico == 0) = 0 ;
% 
% figure;
% box on;
% plot(k,i_inv_estocastico,"r")
% title('Policy Function de la Inversion', tx{:})
% xlabel('Capital $K_t$', tx{:})
% ylabel('Inversion $I_t$', tx{:})
% % legend('Inversion','Location', 'best', tx{:});


%%  2.2 Policy Function
%Creamos vector indicador
vector_indicador = ones(size(Vf_no));
vector_indicador(Vf_inv_estocastico' > Vf_no_inv_estocastico) = 1;
vector_indicador(~(Vf_inv_estocastico' > Vf_no_inv_estocastico)) = 0;

%Ahora creamos inversion
delta = delta1*pdelta1 + delta2*pdelta2;
i_inv_estocastico = (k(kp) - (1-delta)*k);
i_inv_estocastico(vector_indicador == 0) = 0 ;

% Grafique la policy function para esta planta
figure;
box on;
plot(k,i_inv_estocastico,"r")
title('Inversion', tx{:})
xlabel('Capital', tx{:})
title('Policy Function de la Inversion', tx{:})
xlabel('Capital $K_t$', tx{:})
ylabel('Inversion $I_t$', tx{:})

figure;
box on;
% plot(k,(i_inv - i_inv_estocastico),"b")
plot(k,(i_inv_estocastico - i_inv),"b")
title('Inversion', tx{:})
xlabel('Capital', tx{:})
title('Cambios en la Policy Function de la Inversion', tx{:})
xlabel('Capital $K_t$', tx{:})
ylabel('$I_{t, Pregunta\ 2} - I_{t, Pregunta\ 1}$', tx{:})


% Ahora hay que dar la intucion de distintos escenarios. Con alta
% probabilidad que el delta sea alto, eso contribuye a que la inversion sea
% ... Mientras que pasaria lo contrario. %% AHONDAR MAS ACA

%% Pregunta 2.3
% Grafique la funcion Kt+1(Kt) para el escenario en que la depreciacion ha
% sido bajo y para cuando ha sido alta.
% Entonces tenemos que simular 2 veces.
p99 = 0.99;
p01 = 0.01;
delta1 = 0.08;
delta2 = 0.25;
pdelta1 = 0.7;
pdelta2 = 0.3;

[Cti,kti,Vf_inv_estocastico,kp1, iter] = value_inversion_estocastica(beta, alpha, delta1,delta2, A, max_iter, k, tol,pI,F,gamma, p99,p01);
% kp1 ==> posición de capital k+1 óptimo para un kt con inversión con delta aleatorio
[Vf_no_inv_estocastico, policy_choice1, iter1] = value_function_sin_inversion_estocastica(beta, alpha, delta1,delta2,p99,p01, A, max_iter, k, tol);

[Cti,kti,Vf_inv_estocastico,kp2, iter] = value_inversion_estocastica(beta, alpha, delta1,delta2, A, max_iter, k, tol,pI,F,gamma, p01,p99);
% kp2 ==> posición de capital k+1 óptimo para un kt con inversión, con
% delta aleatorio = 25%

[Vf_no_inv_estocastico_2, policy_choice2, iter1] = value_function_sin_inversion_estocastica(beta, alpha, delta1,delta2,p01,p99, A, max_iter, k, tol);

k1 = k*(1-(p99*delta1 + p01*delta2));
% k1 ==> kt+1 con delta = 8%
k2 = k*(1-(p01*delta1 + p99*delta2));
% k2 ==> kt+1 con delta = 25%

% delta = 8%
% pseudo value function con --> (k,k(kp1))
% pseudo value function sin --> (k, k1)

% delta = 25%
% pseudo value function con --> (k,k(kp2))
% pseudo value function sin --> (k, k2)


k(kp1)
figure;
hold on
plot(k, max(k1,k(kp1)),'LineWidth', 4, 'Color', [1 0.5 0.7 1]); % pseudo value function mixta, delta = 8%
plot(k, max(k2,k(kp2)),'LineWidth', 4, 'Color', [1 0.5 0.7 0.5]); % pseudo value function mixta, delta = 25%
plot(k,k(kp1))  % pseudo value function con inversión y delta = 8%
plot(k, k1)
plot(k,k(kp2))  % pseudo value function con inversión y delta = 25%
plot(k, k2)
legend('Probabilidades 1', 'Probabilidades 2', tx{:})
hold off

close all;

figure;
hold on
box on;
plot(k, max(k1,k(kp1)),'LineWidth', 3, 'Color', [0.75 1 0.75 0.75]); % pseudo value function mixta, delta = 8%
plot(k, max(k2,k(kp2)),'LineWidth', 3, 'Color', [0.35 1 0.75 0.5]); % pseudo value function mixta, delta = 25%
plot(k,k(kp1), "b")  % pseudo value function con inversión y delta = 8%
plot(k, k, "r")
plot(k,k(kp2), "m")  % pseudo value function con inversión y delta = 25%
legend(...
    '$K_{t+1}$ con $\delta = 8\%$ y opcion de inversion',...
    '$K_{t+1}$ con $\delta = 25\%$ y opcion de inversion',...
    '$K_{t+1}$ con $\delta = 8\%$ e inversion (sin opcion)',...
    '$K_{t} = K_{t+1}\ \ (45^\circ)$',...
    '$K_{t+1}$ con $\delta = 25\%$ e inversion (sin opcion)',...
    'Location','best',...
    tx{:},...
    'FontSize',10 ...
    )
xlabel('$K_t$', tx{:})
title('Evolucion del capital en funcion del capital vigente', tx{:})
ylabel('$K_{t+1}$', tx{:})
hold off
close all;


%% Pregunta 2.4
% Realice una simulacion de 1000 pasos para esta empresa, comenzando desde
% un capital inicial de K= 0.05. Explique intuituvamente el comportamiento
% ciclico que exhibe esta serie.


T = 1000;  % Number of periods
k_initial = 0.05;  % Initial capital
r = rand(T, 1);
delta = zeros(T, 1);
delta(r < 0.7) = delta1;
delta(r >= 0.7) = delta2;

K = zeros(T, 1);
I = zeros(T, 1);
K(1) = k_initial;


for t = 1:T-1
    [~, idx] = min(abs(k - K(t)));
    i = i_inv_estocastico(idx);
    I(t) = i;
    K(t+1) = (1-delta(t))*K(t) + i;   
end

close all;

% Plot the simulation results
figure;
hold on 
box on;
plot(1:T, K, "Color", "#0072BD");
plot(1:T,I, "Color", "#D14C64");
title('Capital Evolution Over Time', "Interpreter", "Latex");
xlabel('Time Period',"Interpreter", "Latex");
ylabel('Capital Stock',"Interpreter", "Latex");
% legend('Capital', 'Inversion', tx{:},...
%         'Location','best', "Fontsize", 10)
legend('Capital', 'Inversion', tx{:},...
        "Fontsize", 10)
hold off

% Concatenar para desplegar en latex y complementar la interpretación
tabla_informe = [(1:T)' K delta*100 I];
format SHORTG

tabla_informe(1:20,:)
tabla_informe(100:120,:)

format


%% Pregunta 3.1
% Tenemos que obtener la distribuicion ergodica de este problema simulando
% 10.000 veces el problema de la seccion 2.4.  (Guardar ultimo capital e
% inversion)

close all;

T = 1000;  % Number of periods
k_initial = 0.05;  % Initial capital
r = rand(T, 1);
delta = zeros(T, 1);
delta(r < 0.7) = delta1;
delta(r >= 0.7) = delta2;
simulaciones = 10000;
K = zeros(T, 1);
I = zeros(T, 1);
K(1) = k_initial;
distribucion_capital = zeros(1,simulaciones);
distribucion_inversion = zeros(1,simulaciones);
for j=1:simulaciones
    rng(j)
    T = 1000;  % Number of periods
    k_initial = 0.05;  % Initial capital
    r = rand(T, 1);
    delta = zeros(T, 1);
    delta(r < 0.7) = delta1;
    delta(r >= 0.7) = delta2;
    K = zeros(T, 1);
    I = zeros(T, 1);
    for t = 1:T-1
        [~, idx] = min(abs(k - K(t)));
        i = i_inv_estocastico(idx);
        I(t) = i;
        K(t+1) = (1-delta(t))*K(t) + i;   
        
    end
    distribucion_capital(j) = K(end);
    distribucion_inversion(j) = i;
end
% figure;
% hist(distribucion_capital,125);
% figure;
% hist(distribucion_inversion,125);

close all;
figure;

subplot(1,2,1);
box on;
hold on;
histogram(distribucion_capital, 125,'FaceColor',"#8FF8B3",'EdgeColor','none','LineWidth',1.5, 'Normalization','probability')
xlim([0.6 0.95])
title('Distribucion ergodica empirica de $K_T$', "Interpreter", "Latex")
xlabel('$K_T$', "Interpreter", "Latex")
ylabel('Densidad', "Interpreter", "Latex")
hold off;

subplot(1,2,2);
box on;
hold on;
histogram(distribucion_inversion, 125,'FaceColor',"#80B3FF",'EdgeColor','none','LineWidth',1.5, 'Normalization','probability')
xlim([0 0.3])
title('Distribucion ergodica empirica de $I_T$', "Interpreter", "Latex")
xlabel('$I_T$', "Interpreter", "Latex")
ylabel('Densidad', "Interpreter", "Latex")
hold off;


%% Pregunta 3.2


% Intervencion estatal
T = 1000;  % Number of periods
k_initial = 0.05;  % Initial capital
tao = 0.05;
simulaciones = 10000;
K = zeros(T, 1);
I = zeros(T, 1);
K(1) = k_initial;
distribucion_capital = zeros(1,simulaciones);
distribucion_inversion = zeros(1,simulaciones);
for j=1:simulaciones
    rng(j)
    T = 1000;  % Number of periods
    k_initial = 0.05;  % Initial capital
    r = rand(T, 1);
    delta = zeros(T, 1);
    delta(r < 0.7) = delta1 + tao;
    delta(r >= 0.7) = delta2 - tao*(7/3);
    K = zeros(T, 1);
    I = zeros(T, 1);
    for t = 1:T-1
        [~, idx] = min(abs(k - K(t)));
        i = i_inv_estocastico(idx);
        I(t) = i;
        K(t+1) = (1-delta(t))*K(t) + i;   
        
    end
    distribucion_capital(j) = K(end);
    distribucion_inversion(j) = i;
end

close all;
figure;

box on;
histogram(distribucion_capital, 125,'FaceColor',"#8FF8B3",'EdgeColor','none','LineWidth',1.5, 'Normalization','probability')
title('Distribucion ergodica empirica de $K_T$ con impuestos y subsidios', "Interpreter", "Latex")
xlabel('$K_T$', "Interpreter", "Latex")
ylabel('Densidad', "Interpreter", "Latex")
hold off;


delta1 = 0.08+0.05;
delta2 = 0.25-((7/3)*0.05);
pdelta1 = 0.7;
pdelta2 = 0.3;
[Cti,kti,Vf_inv_estocastico,kp, iter] = value_inversion_estocastica(beta, alpha, delta1,delta2, A, max_iter, k, tol,pI,F,gamma, pdelta1,pdelta2);
[Vf_no_inv_estocastico, policy_choice, iter1] = value_function_sin_inversion_estocastica(beta, alpha, delta1,delta2,pdelta1,pdelta2, A, max_iter, k, tol);
figure;
hold on;
box on;
plot(k, Vf_inv_estocastico, "b")
plot(k,Vf_no_inv_estocastico, "r")
Vf_mixta_estocastica = max(Vf_no_inv_estocastico, Vf_inv_estocastico');
plot(k, Vf_mixta_estocastica, 'LineWidth', 8, 'Color', [0.75 1 0.75 0.5]);
title('Value Functions Estocasticas con $\tau$ y $\delta$', tx{:})
xlabel('$Capital\ K_t$', tx{:})
ylabel('$Valor\ Value\ Function$', tx{:})
legend(...
    'V.F. Estocastica Con Inversion $E_0 (V^a(K_t))$',...
    'V.F. Estocastica Sin Inversion $E_0 (V^i(K_t))$',...
    'V.F. Estocastica Con Opcion de Inversion $E_0(V(K_t))$',...
    'Location','best',...
    tx{:},...
    'FontSize',10 ...
        );
hold off;

%%  3.2 Policy Function
%Creamos vector indicador
vector_indicador = ones(size(Vf_no));
vector_indicador(Vf_inv_estocastico' > Vf_no_inv_estocastico) = 1;
vector_indicador(~(Vf_inv_estocastico' > Vf_no_inv_estocastico)) = 0;

%Ahora creamos inversion
delta = delta1*pdelta1 + delta2*pdelta2;
i_inv_estocastico = (k(kp) - (1-delta)*k);
i_inv_estocastico(vector_indicador == 0) = 0 ;

% Grafique la policy function para esta planta
figure;
box on;
plot(k,i_inv_estocastico,"r")
title('Inversion', tx{:})
xlabel('Capital', tx{:})
title('Policy Function de la Inversion con $\tau$ y $\delta$', tx{:})
xlabel('Capital $K_t$', tx{:})
ylabel('Inversion $I_t$', tx{:})

close all;

%% Rincon de Funciones 
function [prod] = produccion(capital,alpha)
    prod = capital.^alpha;
end

% Funcion de iteracion P1 
function [Cti,kti,Vf,kp, iter] = value_inversion(beta, alpha, delta, A, max_iter, k, tol,pI,F,gamma)
e = 1; % Error inicial
iter = 0;
Vf = ones(length(k),1); %Guess inicial

while e > tol && iter < max_iter
    %Prealocation
    Vp = NaN(length(k),1);  
    kp = NaN(length(k),1);  
    Caux = A*k'.^alpha - pI*(k - (1-delta)*k') - F*k' - (gamma/2).*((k - (1-delta)*k')./k').^2.*k';
    % k es un vector fila, al restar, los valor de k se restan en filas a
    % varias columnas pegas de k'. 
    Caux(Caux<=0) = NaN;  
    Vaux= Caux + beta*Vf';    % ecuación de bellman
    [Vp(:,1),kp(:,1)]= max(Vaux,[],2);                           
    e = max(abs(Vf-Vp));    
    iter = iter + 1 ;
    Vf = Vp;
end

Cti = A.* k.^alpha + (1- delta).*k - k(kp');%Consumption policy
kti = k(kp');%Assets policy
Vf = Vf';
kp = kp';
end



function [Vf, policy_choice, iter] = value_function_sin_inversion(beta, alpha, delta, A, max_iter, k_grid, tol)

    % Inicializacion
    e = 1; % Error Inicial
    iter = 0; % Iteraciones
    Vf = zeros(length(k_grid), 1); % Guess Inicial
    policy_choice = zeros(length(k_grid), 1); % Eleccion Optima

    % Make sure k_grid is a column vector
    k_grid = k_grid(:);

    % Precompute repeated values
    k_grid_depreciated = k_grid * (1 - delta);

    while e > tol && iter < max_iter
        % Computar funcion de valor de interpolar
        Vi = A * k_grid.^alpha + beta * interp1(k_grid, Vf, k_grid_depreciated, 'linear', 'extrap');
        % Elegir el maximo
        [Vp, policy_indices] = max(Vi, [], 2);
        
        % Computamos error y agregamos iteracion
        e = max(abs(Vf - Vp));
        iter = iter + 1;
        Vf = Vp; % Actualizar value function
        policy_choice = policy_indices; % Actualizar policy
    end

end

% Funcion de iteracion P1 
function [Cti,kti,Vf,kp, iter] = value_inversion_estocastica(beta, alpha, delta1,delta2, A, max_iter, k, tol,pI,F,gamma, pdelta1,pdelta2)
e = 1; % Error inicial
iter = 0;
Vf = ones(length(k),1); %Guess inicial

while e > tol && iter < max_iter
    %Prealocation
    Vp = NaN(length(k),1);  
    kp = NaN(length(k),1);  
    Caux1 = A*k'.^alpha - pI*(k - (1-delta1)*k') - F*k' - (gamma/2).*((k - (1-delta1)*k')./k').^2.*k';
    Caux2 = A*k'.^alpha - pI*(k - (1-delta2)*k') - F*k' - (gamma/2).*((k - (1-delta2)*k')./k').^2.*k';
    Caux = pdelta1*Caux1 +pdelta2*Caux2;
    % k es un vector fila, al restar, los valor de k se restan en filas a
    % varias columnas pegas de k'. 
    Caux(Caux<=0) = NaN;  
    Vaux= Caux + beta*Vf';    % ecuación de bellman
    [Vp(:,1),kp(:,1)]= max(Vaux,[],2);     
    e = max(abs(Vf-Vp));                          
    iter = iter + 1 ;
    Vf = Vp;
end

Cti = A.* k.^alpha + (1- delta1).*k - k(kp');%Consumption policy
kti = k(kp');%Assets policy
Vf = Vf';
kp = kp';
end

function [Vf, policy_choice, iter] = value_function_sin_inversion_estocastica(beta, alpha, delta1,delta2,pdelta1,pdelta2, A, max_iter, k_grid, tol)

    % Inicializacion
    e = 1; % Error Inicial
    iter = 0; % Iteraciones
    Vf = zeros(length(k_grid), 1); % Guess Inicial
    policy_choice = zeros(length(k_grid), 1); % Eleccion Optima

    % Make sure k_grid is a column vector
    k_grid = k_grid(:);

    % Precompute repeated values
    k_grid_depreciated = (k_grid * (1 - delta1))*pdelta1 + (k_grid * (1 - delta2))*pdelta2;

    while e > tol && iter < max_iter
        % Computar funcion de valor de interpolar
        Vi  = A * k_grid.^alpha + beta * interp1(k_grid, Vf, k_grid_depreciated, 'linear', 'extrap');
        % Elegir el maximo
        [Vp, policy_indices] = max(Vi, [], 2);
        
        % Computamos error y agregamos iteracion
        e = max(abs(Vf - Vp));
        iter = iter + 1;
        Vf = Vp; % Actualizar value function
        policy_choice = policy_indices; % Actualizar policy
    end

end