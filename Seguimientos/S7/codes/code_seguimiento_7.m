%% Seguimiento 7

clear; 
close all;
clc;
seed = 1404;
rng(seed)

% Parámetros
beta = 0.9615;
R = 1.04; % R = (1 + r) = (1 + 0.04)
beta = 0.9615;
R = 1.04; % R = (1 + r) = (1 + 0.04)
sigma = 2;
phi = 1;
w = 8;
Amin = -phi;
Amax = 30;
N_A = 100;
gridA = linspace(Amin,Amax,N_A);
N_Z = 2;
Z_l = 0.5709;
Z_h = 1.4291;
gridZ = linspace(Z_l, Z_h, N_Z);
%shocks
theta = w.*gridZ'; % Senda de ingresos para el estado alto y bajo
k = 0.9371;
P = [k 1-k; 1-k k];


%% Pregunta 1
[CC]= endo_infinite (gridA, N_A, gridZ, N_Z, w, P, beta, sigma, R, phi);

policy = CC(:,:,size(CC, 3));
policy_bad_state = policy(1,:);
policy_good_state = policy(2,:);

figure;
plot(gridA,policy_good_state,"g");
hold on;
box on;
title('Converged consumption', 'Interpreter','Latex','FontSize', 15)
plot(gridA,policy_bad_state,"r");
xlabel('$Activos\ A_t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Consumo\ C_t$', 'Interpreter','Latex','FontSize', 15)
legend(...
    'Policy cuando el salario es $w_t = 11.4328$',...
    'Policy cuando el salario es $w_t = 4.5672$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([Amin Amax])
hold off;
saveas(gcf,"policy.svg")
close all;


%% Pregunta 2

T = 1000;
estado_actual=round(rand(1))+1;
ahorro_inicial = 10;

display(estado_actual)  % == 2
historial_estados = Markov(P, estado_actual, T)';

historial_consumo = NaN(1,T);
historial_ahorro = NaN(1,T);
historial_ingresos = NaN(1,T);

historial_ahorro(1) = ahorro_inicial;
ahorro_actual = ahorro_inicial;

for i = 1:T

    historial_consumo(i) = interp1(gridA, CC(historial_estados(i),:,size(CC, 3)), ahorro_actual,[],'extrap');
    historial_ingresos(i) = theta(historial_estados(i));
    if i == T
        continue;
    else
        historial_ahorro(i+1) = (R)*ahorro_actual + historial_ingresos(i) - historial_consumo(i);
        % Restricción de deuda
        if historial_ahorro(i+1) < -phi
            historial_ahorro(i+1) = -phi; 
            historial_consumo(i) = (R)*ahorro_actual + historial_ingresos(i) - historial_ahorro(i+1);
        end
        % Guardamos este ahorro de mañana para iterar nuevamente para
        % el dia siguiente
        ahorro_actual = historial_ahorro(i+1);
    end
end

plot(1:T, historial_consumo)
hold on;
plot(1:T, historial_ahorro)
plot(1:T, historial_ingresos)
hold off;


figure;
plot(1:T,historial_consumo,"g");
hold on;
box on;
plot(1:T,historial_ahorro,"b");
plot(1:T,historial_ingresos,"k");
title('Simulation', 'Interpreter','Latex','FontSize', 15)
xlabel('$t$', 'Interpreter','Latex','FontSize', 15)
ylabel('$Realizaciones$', 'Interpreter','Latex','FontSize', 15)
legend(...
    '$Consumo$',...
    '$Ahorro$',...
    '$Ingresos$',...
    'Location','best',...
    'Interpreter','Latex',...
    'FontSize',10 ...
        );
xlim([1 T])
hold off;
saveas(gcf,"simulation.svg")
close all;


