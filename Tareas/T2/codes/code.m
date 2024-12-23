%% Tarea 1

clear; 
close all;
clc;
seed = 1404;
rng(seed)


% Parámetros
beta = 0.95;
tol = 0.0001;
Kmin = 0.05;
Kmax = 200;
Nk = 1000;
A = 1.25;
alpha = 0.28;
gamma = 1/3;
Delta = 0.08;
a = Kmax;


gridK = linspace(Kmin,Kmax,Nk);


%% Pregunta 1

[II]= endo_infinite (gridK, beta, alpha, a, A, Delta, gamma, tol, Kmin, Kmax);

tx= {'Interpreter','Latex','FontSize', 15};

figure;
hold on;
box on 
plot(gridK,II(:,size(II, 2)), "g")
title('Policy Function', tx{:})
xlabel('$Capital\ K_t$', tx{:})
ylabel('$Inversion\ I_t$', tx{:})
hold off;
saveas(gcf,"1_c.1_policy_inversion.svg")


% Policy del capital

K_t_pus_1 = (1-Delta)*gridK + II(:,size(II, 2))';

figure;
hold on;
box on 
plot(gridK,K_t_pus_1, "b")
title('Policy Function', tx{:})
xlabel('$Capital\ K_t$', tx{:})
ylabel('$Capital\ K_{t+1}$', tx{:})
hold off;
saveas(gcf,"1_c.2_policy_capital.svg")


% Simulación para K0 = 0.05

T = 100;

historial_capital = zeros(T,1);
historial_inversion = zeros(T,1);
historial_capital(1) = 0.05;
historial_inversion(1) = interp1(gridK,II(:,size(II, 2)),historial_capital(1),[],'extrap');
for t = 2:T
    historial_capital(t) = (1-Delta)*historial_capital(t-1) + historial_inversion(t-1);
    historial_inversion(t) = interp1(gridK,II(:,size(II, 2)),historial_capital(t),[],'extrap');
end


figure;
plot(1:T,historial_capital,"b");
hold on;
box on;
plot(1:T,historial_inversion,"g");
title('Simulation $K_0 = 0.05$', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones Simulacion$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Capital$',...
    '$Inversion$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"1_d.1_simulation_K_0.05.svg")
close all;

kap_estacionario_K_0_05 = historial_capital(T)
inv_estacionaria_K_0_05 = historial_inversion(T)


% Simulación para K0 = 90


historial_capital = zeros(T,1);
historial_inversion = zeros(T,1);
historial_capital(1) = 90;
historial_inversion(1) = interp1(gridK,II(:,size(II, 2)),historial_capital(1),[],'extrap');
for t = 2:T
    historial_capital(t) = (1-Delta)*historial_capital(t-1) + historial_inversion(t-1);
    historial_inversion(t) = interp1(gridK,II(:,size(II, 2)),historial_capital(t),[],'extrap');
end


figure;
plot(1:T,historial_capital,"b");
hold on;
box on;
plot(1:T,historial_inversion,"g");
title('Simulation $K_0 = 90$', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones Simulacion$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Capital$',...
    '$Inversion$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"1_d.2_simulation_K_90.svg")
close all;

kap_estacionario_K_90 = historial_capital(T)
inv_estacionaria_K_90 = historial_inversion(T)


% Simulación para K0 = 190


historial_capital = zeros(T,1);
historial_inversion = zeros(T,1);
historial_capital(1) = 190;
historial_inversion(1) = interp1(gridK,II(:,size(II, 2)),historial_capital(1),[],'extrap');
for t = 2:T
    historial_capital(t) = (1-Delta)*historial_capital(t-1) + historial_inversion(t-1);
    historial_inversion(t) = interp1(gridK,II(:,size(II, 2)),historial_capital(t),[],'extrap');
end


figure;
plot(1:T,historial_capital,"b");
hold on;
box on;
plot(1:T,historial_inversion,"g");
title('Simulation $K_0 = 190$', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones Simulacion$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Capital$',...
    '$Inversion$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"1_d.3_simulation_K_190.svg")
close all;

kap_estacionario_K_190 = historial_capital(T)
inv_estacionaria_K_190 = historial_inversion(T)


% Capitales e Inversiones Estacionarias


tabla_aux_para_latex = [...
    [0.05 90 190]'...
    [kap_estacionario_K_0_05 kap_estacionario_K_90 kap_estacionario_K_190]'...
    [inv_estacionaria_K_0_05 inv_estacionaria_K_90 inv_estacionaria_K_190]'...
    ]


%% Pregunta 2

save('workspace_a_la_2')
clear all;
load('workspace_a_la_2')

[II_cournot]= endo_infinite_cournot (gridK, beta, alpha, a, A, Delta, gamma, tol, Kmin, Kmax);

grid_cournot_k_i = repmat(gridK', 1, size(gridK, 2));
grid_cournot_k_j = repmat(gridK, size(gridK, 2), 1);

save('workspace_a_la_2_c')
clear all;
load('workspace_a_la_2_c')

% format SHORTG


% set(gcf,'Renderer','Painters')
plot_policy = mesh(grid_cournot_k_i,grid_cournot_k_j,II_cournot(:,:,size(II_cournot, 3)),'FaceAlpha','0.5')
plot_policy.FaceColor = 'flat';
title("Policy Function en funcion de capitales $K_{t+1}^i, K_{t+1}^j$",'interpreter','latex')
xlabel('$K_{t+1}^i$','interpreter','latex')
ylabel('$K_{t+1}^j$','interpreter','latex')
zlabel('$I_t^i$','interpreter','latex')
saveas(gcf,"1.2_0_policy_3d.svg")


%% Pregunta 1.2.c

matriz_de_convergencia = II_cournot(:,:,size(II_cournot, 3));

% j chica

policy_contra_chica = matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<0.05)+1);

% format SHORTG

figure;
hold on;
box on 
plot(gridK,policy_contra_chica, "g")
title('Policy Function contra firma chica', tx{:})
xlabel('$Capital\ K_t^i$', tx{:})
ylabel('$Inversion\ I_t$', tx{:})
hold off;
saveas(gcf,"1.2_c.1_policy_inversion_contra_chica.svg")


% j mediana

policy_contra_mediana = matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<90)+1);

figure;
hold on;
box on 
plot(gridK,policy_contra_mediana, "g")
title('Policy Function contra firma mediana', tx{:})
xlabel('$Capital\ K_t^i$', tx{:})
ylabel('$Inversion\ I_t$', tx{:})
hold off;
saveas(gcf,"1.2_c.2_policy_inversion_contra_mediana.svg")


% j grande

policy_contra_grande = matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<190)+1);

figure;
hold on;
box on 
plot(gridK,policy_contra_grande, "g")
title('Policy Function contra firma grande', tx{:})
xlabel('$Capital\ K_t^i$', tx{:})
ylabel('$Inversion\ I_t$', tx{:})
hold off;
saveas(gcf,"1.2_c.3_policy_inversion_contra_grande.svg")


% Policies juntas

close all;
figure;
hold on;
box on 
plot(gridK,policy_contra_chica, "g")
plot(gridK,policy_contra_mediana, "b")
plot(gridK,policy_contra_grande, "r")
legend(...
    '$Policy\ Function\ contra\ firma\ chica$',...
    '$Policy\ Function\ contra\ firma\ mediana$',...
    '$Policy\ Function\ contra\ firma\ grande$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
title('Policy Functions', tx{:})
xlabel('$Capital\ K_t^i$', tx{:})
ylabel('$Inversion\ I_t$', tx{:})
hold off;
saveas(gcf,"1.2_c.4_policy_functions_juntas.svg")

% [policy_contra_chica policy_contra_mediana policy_contra_grande]

% format


K_t_pus_1_chica = (1-Delta)*gridK' + policy_contra_chica
K_t_pus_1_mediana = (1-Delta)*gridK' + policy_contra_mediana
K_t_pus_1_grande = (1-Delta)*gridK' + policy_contra_grande

figure;
hold on;
box on 
plot(gridK,K_t_pus_1_chica, "g")
plot(gridK,K_t_pus_1_mediana, "b")
plot(gridK,K_t_pus_1_grande, "r")
legend(...
    '$Policy\ Function\ contra\ firma\ chica$',...
    '$Policy\ Function\ contra\ firma\ mediana$',...
    '$Policy\ Function\ contra\ firma\ grande$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
title('Policy Function', tx{:})
xlabel('$Capital\ K_t$', tx{:})
ylabel('$Capital\ K_{t+1}$', tx{:})
hold off;
saveas(gcf,"1.2_c.5_policy_capital_juntas.svg")




%% Pregunta 1.2.d

[kap_estacionario_K_0_05 kap_estacionario_K_90 kap_estacionario_K_190]


% Simulación para K_j_0 = 0.05

T = 100;

historial_capital_i = zeros(T,1);
historial_inversion_i = zeros(T,1);

historial_capital_j = zeros(T,1);
historial_inversion_j = zeros(T,1);

historial_capital_i(1) = kap_estacionario_K_0_05;
historial_inversion_i(1) = interp1(gridK,policy_contra_chica,historial_capital_i(1),[],'extrap');


historial_capital_j(1) = 0.05;
%%%%%%%%%%%%%%%

historial_inversion_j(1) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<kap_estacionario_K_0_05)+1),...
    historial_capital_j(1),...
    [],'extrap');


for t = 2:T
    historial_capital_i(t) = (1-Delta)*historial_capital_i(t-1) + historial_inversion_i(t-1);
    historial_capital_j(t) = (1-Delta)*historial_capital_j(t-1) + historial_inversion_j(t-1);
    
    % Policy i
    historial_inversion_i(t) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<historial_capital_j(t))+1),...
        historial_capital_i(t),...
        [],'extrap');
    
    % Policy j
    historial_inversion_j(t) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<historial_capital_i(t))+1),...
        historial_capital_j(t),...
        [],'extrap');
    
end

% Simulación para K_j_0 = 90

historial_capital_i_90 = zeros(T,1);
historial_inversion_i_90 = zeros(T,1);

historial_capital_j_90 = zeros(T,1);
historial_inversion_j_90 = zeros(T,1);

historial_capital_i_90(1) = kap_estacionario_K_90;
historial_inversion_i_90(1) = interp1(gridK,policy_contra_mediana,historial_capital_i_90(1),[],'extrap');


historial_capital_j_90(1) = 90;
%%%%%%%%%%%%%%%

historial_inversion_j_90(1) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<kap_estacionario_K_90)+1),...
    historial_capital_j_90(1),...
    [],'extrap');


for t = 2:T
    historial_capital_i_90(t) = (1-Delta)*historial_capital_i_90(t-1) + historial_inversion_i_90(t-1);
    historial_capital_j_90(t) = (1-Delta)*historial_capital_j_90(t-1) + historial_inversion_j_90(t-1);
    
    % Policy i
    historial_inversion_i_90(t) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<historial_capital_j_90(t))+1),...
        historial_capital_i_90(t),...
        [],'extrap');
    
    % Policy j
    historial_inversion_j_90(t) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<historial_capital_i_90(t))+1),...
        historial_capital_j_90(t),...
        [],'extrap');
    
end

% Simulación para K_j_0 = 190

historial_capital_i_190 = zeros(T,1);
historial_inversion_i_190 = zeros(T,1);

historial_capital_j_190 = zeros(T,1);
historial_inversion_j_190 = zeros(T,1);

historial_capital_i_190(1) = kap_estacionario_K_190;
historial_inversion_i_190(1) = interp1(gridK,policy_contra_mediana,historial_capital_i_190(1),[],'extrap');


historial_capital_j_190(1) = 190;
%%%%%%%%%%%%%%%

historial_inversion_j_190(1) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<kap_estacionario_K_190)+1),...
    historial_capital_j_190(1),...
    [],'extrap');


for t = 2:T
    historial_capital_i_190(t) = (1-Delta)*historial_capital_i_190(t-1) + historial_inversion_i_190(t-1);
    historial_capital_j_190(t) = (1-Delta)*historial_capital_j_190(t-1) + historial_inversion_j_190(t-1);
    
    % Policy i
    historial_inversion_i_190(t) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<historial_capital_j_190(t))+1),...
        historial_capital_i_190(t),...
        [],'extrap');
    
    % Policy j
    historial_inversion_j_190(t) = interp1(gridK,...
    matriz_de_convergencia(:,sum(grid_cournot_k_j(1,:)<historial_capital_i_190(t))+1),...
        historial_capital_j_190(t),...
        [],'extrap');
    
end


figure;
plot(1:T,historial_capital_i,"g");
hold on;
plot(1:T,historial_capital_i_90,"b");
plot(1:T,historial_capital_i_190,"r");
box on;
title('Simulations vs $j$', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones\ Capital\ i$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Capital_j = 0.05$',...
    '$Capital_j = 90$',...
    '$Capital_j = 190$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"1.2_d.1_simulation_K_de_i.svg")
close all;

figure;
plot(1:T,historial_inversion_i,"g");
hold on;
plot(1:T,historial_inversion_i_90,"b");
plot(1:T,historial_inversion_i_190,"r");
box on;
title('Simulations vs $j$', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones\ Inversion\ i$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Inversion\ Capital_j = 0.05$',...
    '$Inversion\ Capital_j = 90$',...
    '$Inversion\ Capital_j = 190$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"1.2_d.2_simulation_I_de_i.svg")
close all;

figure;
plot(1:T,historial_capital_j,"g");
hold on;
plot(1:T,historial_capital_j_90,"b");
plot(1:T,historial_capital_j_190,"r");
box on;
title('Simulations vs $i$', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones\ Capital\ j$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Capital_j = 0.05$',...
    '$Capital_j = 90$',...
    '$Capital_j = 190$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"1.2_d.3_simulation_K_de_j.svg")
close all;

figure;
plot(1:T,historial_inversion_j,"g");
hold on;
plot(1:T,historial_inversion_j_90,"b");
plot(1:T,historial_inversion_j_190,"r");
box on;
title('Simulations vs $i$', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones\ Inversion\ j$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Inversion\ Capital_j = 0.05$',...
    '$Inversion\ Capital_j = 90$',...
    '$Inversion\ Capital_j = 190$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"1.2_d.3_simulation_I_de_j.svg")
close all;



kap_estacionario_Kj_0_05_cournot = historial_capital_i(T)
inv_estacionaria_Kj_0_05_cournot = historial_inversion_i(T)

kap_estacionario_Kj_90_cournot = historial_capital_i_90(T)
inv_estacionaria_Kj_90_cournot = historial_inversion_i_90(T)

kap_estacionario_Kj_190_cournot = historial_capital_i_190(T)
inv_estacionaria_Kj_190_cournot = historial_inversion_i_190(T)




%% Pregunta 1.2 e

% El resultado final debería ser algo como esto

tabla_aux_2_para_latex = [...
    [0.05 90 190]'...
    [kap_estacionario_K_0_05 kap_estacionario_K_90 kap_estacionario_K_190]'...
    [kap_estacionario_Kj_0_05_cournot kap_estacionario_Kj_90_cournot kap_estacionario_Kj_190_cournot ]'...
    [inv_estacionaria_K_0_05 inv_estacionaria_K_90 inv_estacionaria_K_190]'...
    [inv_estacionaria_Kj_0_05_cournot inv_estacionaria_Kj_90_cournot inv_estacionaria_Kj_190_cournot]'...
    ]


