%% Pregunta 1.1

clear all; close all;

% Parametros
alpha = 0.5;
M = 100;
p1 = 15;
p2 = 10;
max_x1 = M/p1;
max_x2 = M/p2;
N_grid = 100;

% Grilla de consumos posibles
x1_grid = linspace(0,max_x1,N_grid)';  % linspace(X1, X2, N): N points between X1 and X2
x2_grid = linspace(0,max_x2,N_grid)';

% 10 primeras y últimas observaciones
display_consumos = [(1:length(x1_grid))' x1_grid x2_grid];
% id, c1, c2
[display_consumos(1:10,:); display_consumos(length(display_consumos)-9:end,:)]


% Replicación manual de meshgrid para la Cobb Douglas
X1 = repmat(linspace(0,max_x1,N_grid), N_grid, 1);
X2 = repmat(linspace(0,max_x2,N_grid)', 1, N_grid);

utilidad = udd_cobb_douglas(alpha, X1, X2);

set(gcf,'Renderer','Painters')
plot_f_de_u = mesh(X1,X2,utilidad,'FaceAlpha','0.5')
plot_f_de_u.FaceColor = 'flat';
title("Funci\'on de Utilidad $U(x_1, x_2) = {x_1}^{0.5}{x_2}^{0.5}$",'interpreter','latex')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$utils$','interpreter','latex')



%% Pregunta 1.2


gasto = funcion_de_gasto(p1,p2,X1,X2);
figure
set(gcf,'Renderer','Painters')
plot_f_de_gasto = mesh(X1,X2,gasto,'FaceAlpha','0.5')
plot_f_de_gasto.FaceColor = 'flat';
title("Funci\'on de Gasto $G(p_1, p_2, x_1, x_2) = p_1 x_1+ p_2 x_2$",'interpreter','latex')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$G$','interpreter','latex')


[zonas_de_gasto, parametros_plot_gasto] = contourf(...
    X1,X2,gasto,[0:20:200],...
     'ShowText','on'...
    )
colorbar
set(parametros_plot_gasto,'LineColor','flat')
clabel(zonas_de_gasto,parametros_plot_gasto,'FontSize',15,'Color','w','FontName','Courier')
title("Zonas de Gasto = $\{(x_1, x_2): x_1 \geq 0, x_2 \geq 0, p_1 x_1+ p_2 x_2 \leq Level\}$",'interpreter','latex')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')


%% Pregunta 1.3

% Grilla de consumos posibles con RP con igualdad...

% x1_grid corresponde a los posibles consumos de x1;
% El agente tiene una RP vector que va desde todo el presupuesto
% gastado en x1, a todo el presupuesto en x2
% % p1*x1 + p2*x2 = M ==> x2 = (1/p2)*(M - p1*x1)
x2_grid = (1/p2)*(M - p1*x1_grid);

%  Cálculo de la máxima utilidad...

utilidad_RP = NaN(length(x1_grid), 1);
for i=1:length(x1_grid)
     utilidad_RP(i) = udd_cobb_douglas_para_escalares(alpha, x1_grid(i), x2_grid(i));  
end
max_utilidad = max(utilidad_RP);

% Cálculo del punto tangete entre curva de
% indiferencia en máx. udd y la RP
optimal_consumption_index = find(utilidad_RP==max_utilidad);

% Hay dos puntos óptimos

% Punto tangente
A_x1 = x1_grid(optimal_consumption_index(1))
A_x2 = x2_grid(optimal_consumption_index(1))

% Sólo para ver el segundo punto óptimo
x1_grid(optimal_consumption_index(2))
x2_grid(optimal_consumption_index(2))

% Curvas de nivel de la función de utilidad

% Gráfico con contourf

[curvas_indiferencia, parametros_plot_indiferencia] = contourf(...
    X1,X2,utilidad,[linspace(0,2*max_utilidad,11)],'ShowText','on'...
)
colorbar
set(parametros_plot_indiferencia,'LineColor','flat')
clabel(curvas_indiferencia,parametros_plot_indiferencia,...
    'FontSize',15,'Color','w','FontName','Courier')
title("Curvas de indiferencia Cobb-Douglas",'interpreter','latex')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
hold on;
plot(...
    A_x1,A_x2,...
    'o-','MarkerFaceColor','green','MarkerEdgeColor','green',...
    'MarkerSize',6 ...
)
text(A_x1+.25,A_x2+.25,'A','interpreter','latex','Color','green','FontSize',15)
plot(x1_grid,x2_grid,'g','LineWidth',1)
hold off;


% Gráfico con contour

plot(x1_grid, x2_grid, 'r')
hold on;
[curvas_indiferencia, parametros_plot_indiferencia] = contour(...
    X1,X2,utilidad,[linspace(0,2*max_utilidad,11)],'ShowText','on'...
)
colorbar
set(parametros_plot_indiferencia,'LineColor','flat')
clabel(curvas_indiferencia,parametros_plot_indiferencia,...
    'FontSize',15,'Color','b','FontName','Courier')
title("Curvas de indiferencia Cobb-Douglas",'interpreter','latex')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
hold on;
plot(...
    A_x1,A_x2,...
    'o-','MarkerFaceColor','red','MarkerEdgeColor','red',...
    'MarkerSize',6 ...
)
text(A_x1+.25,A_x2+.25,'A','interpreter','latex','Color','red','FontSize',14)
hold off;



%% Pregunta 2.1.1

clear all; close all;

% Parametros
A = 1;
alpha = 0.25;
r = 1.05;
p = 1;
L = 2;
N_grid = 100;
max_k = 1;

% Grilla de capitales posibles
k_grid = linspace(0,max_k,N_grid)';  % linspace(X1, X2, N): N points between X1 and X2

profits = beneficios(p,A,k_grid,alpha,L,r);
maximo_beneficio = max(profits);

% Beneficios de las 10 primeras y últimas observaciones de capital discreto
display_beneficios = [(1:length(k_grid))' k_grid profits];
[display_beneficios(1:10,:); display_beneficios(length(display_beneficios)-9:end,:)]


capital_optimo = k_grid(find(profits == maximo_beneficio));

capital_optimo
maximo_beneficio

% Plot
plot(k_grid, profits, 'Color', 'black');
set(gca,'xticklabel',num2str(''),'yticklabel',num2str(''))
hold on; 

ylab = ylabel("$$\Pi$$", 'FontSize', 15, 'Position',[0 0.8], 'Interpreter', 'latex');
set(ylab,'rotation',0,'VerticalAlignment','bottom')

xlabel("$$K_t$$", 'FontSize', 15, 'Position', [1 -0.01], 'Interpreter', 'latex') 
set(gca, 'Box', 'off');

plot(capital_optimo, maximo_beneficio,'.','MarkerSize', 20, 'Color', 'black')
line(...
    [capital_optimo capital_optimo],...
    [0 maximo_beneficio], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 0.6...
);
line(...
    [0 capital_optimo],...
    [maximo_beneficio maximo_beneficio], 'LineStyle', '--', 'Color', 'black', 'LineWidth', 0.6...
);
text(capital_optimo - 0.125, - 0.025, '$$K_t^*=0.2323$$','FontSize', 9, 'Interpreter', 'latex')
text( - 0.16, maximo_beneficio - 0.00125, '$$\Pi^*=0.7379$$','FontSize', 9, 'Interpreter', 'latex')
hold off;


%% Pregunta 2.1.2

min_r = 1;
max_r = 1.2;
N_r_grid = 50;

% Dda Kt = (r/(\alpha * p * A * \overline{L}^0.5))^(1/\alpha - 1)

r_grid = linspace(min_r,max_r,N_r_grid)';
Kt = (r_grid*(1/(alpha * p * A * L^0.5))).^(1/(alpha - 1));

% Plot
plot(r_grid, Kt, 'Color', 'black');
set(gca, 'Box', 'off');
hold on;

ylab = ylabel("$$r$$", 'FontSize', 15, 'Position',[1 0.2515], 'Interpreter', 'latex');
set(ylab,'rotation',0,'VerticalAlignment','bottom')

xlabel("$$K_t$$", 'FontSize', 15, 'Position', [1.21 0.192], 'Interpreter', 'latex') 
hold off;


%% Pregunta 2.2.1

close all;
rng(73);  % Seed
N_sim  = 10000;

% Parámetros tabla enunciado
beta = 1;
r = 1.05;
z1 = 10;
rho = 0.88;
N_b_grid = 1000;

epsilon = normrnd(0,sqrt(1/2),N_sim,1);
z2 = rho * ones(length(epsilon), 1) * sqrt(z1) + epsilon;

% 10 primeras y últimas observaciones de la simulación
display_z2 = [(1:length(z2))' z2];
format SHORTG
[display_z2(1:10,:); display_z2(length(display_z2)-9:end,:)]
format short


b_optimo_por_sim = NaN(N_sim, 1);

for i=1:N_sim
    b_grid = linspace(-z2(i)/r, z1, N_b_grid+2)';
    b_grid = b_grid(2:end-1);
    max_b = find(...
        log(z1 - b_grid) + beta * log(b_grid*r + z2(i)) == max(log(z1 - b_grid) + beta * log(b_grid*r + z2(i)))...
        );
    if length(max_b) ~= 1  % Si hay más de un máximo
        b_optimo_por_sim(i) = b_grid(max_b(1));  % toma el primero (supuesto)
    else
        b_optimo_por_sim(i) = b_grid(max_b);
    end
end

mean(b_optimo_por_sim)

histogram_bs = histfit(b_optimo_por_sim);
histogram_bs(1).EdgeColor = "none";
histogram_bs(1).FaceColor = "#A1FEC1";
histogram_bs(2).Color = "#00F0EC";
histogram_bs(2).LineWidth = 0.75;
xlabel('$b^*$', 'Interpreter', 'latex');
xlim([2.5 5])
ylab = ylabel("$Frecuencia$", 'Interpreter', 'latex');
hold off;


%% Pregunta 2.2.2

b_promedio_optimo_por_r_movil = NaN(N_r_grid, 1);


for r_movil=linspace(1, 1.2, N_r_grid);
    
    % Función de Montecarlo
    z2 = montecarlo_ingresos(0,sqrt(1/2),N_sim,rho,z1);
    
    for i=1:N_sim
        b_grid = linspace(-z2(i)/r_movil, z1, N_b_grid+2)';
        b_grid = b_grid(2:end-1);
        max_b = find(...
            log(z1 - b_grid) + beta * log(b_grid*r_movil + z2(i)) == max(log(z1 - b_grid) + beta * log(b_grid*r_movil + z2(i)))...
            );
        if length(max_b) ~= 1  % Si hay más de un máximo
            b_optimo_por_sim(i) = b_grid(max_b(1));  % toma el primero (supuesto)
        else
            b_optimo_por_sim(i) = b_grid(max_b);
        end
    end
    
    b_promedio_optimo_por_r_movil(sum(isnan(b_promedio_optimo_por_r_movil)==0)+1) = mean(b_optimo_por_sim);

end

b_promedio_optimo_por_r_movil;

% Plot
plot(r_grid, b_promedio_optimo_por_r_movil, 'Color', 'black');
set(gca, 'Box', 'off');
hold on; 

ylab = ylabel("$$r$$", 'FontSize', 15, 'Position',[1 3.855], 'Interpreter', 'latex');
set(ylab,'rotation',0,'VerticalAlignment','bottom')

xlabel("$$b^*$$", 'FontSize', 15, 'Position', [1.21 3.61], 'Interpreter', 'latex') 
hold off;


%% 2.2.3

b_promedio_optimo_por_r_movil_fsolve = NaN(N_r_grid, 1);

count_of_r_s = 1;
for r_movil=linspace(1, 1.2, N_r_grid);
    % Para cada r

    % Función de Montecarlo
    z2 = montecarlo_ingresos(0,sqrt(1/2),N_sim,rho,z1);
    z2 = mean(z2);

    umg = @(b_star) ((-1/(z1 - b_star)) + beta*r_movil*(1/(b_star*r_movil + z2)))
    b_star = fsolve(umg, 3.5)
    b_promedio_optimo_por_r_movil_fsolve(count_of_r_s) = b_star

    count_of_r_s = count_of_r_s + 1;
    
end
% Checkpoint

% Plot
plot(r_grid, b_promedio_optimo_por_r_movil_fsolve, 'Color', 'black');
set(gca, 'Box', 'off');
hold on; 

ylab = ylabel("$$r$$", 'FontSize', 15, 'Position',[1 3.855], 'Interpreter', 'latex');
set(ylab,'rotation',0,'VerticalAlignment','bottom')

% Cambiar position
xlabel("$$b^*$$", 'FontSize', 15, 'Position', [1.21 3.61], 'Interpreter', 'latex') 
hold off;


%% 2.2.4

% Parametros
alpha = 0.45;
z1 = 10;
rho = 0.88;
beta = 0.85;
p = 1.15;
A = 3;
L = 2;
N_r_grid = 50;

% Demanda

min_r = 1;
max_r = 1.2;
N_r_grid = 50;

% Dda Kt = (r/(\alpha * p * A * \overline{L}^0.5))^(1/\alpha - 1)

r_grid = linspace(min_r,max_r,N_r_grid)';
Kt = (r_grid*(1/(alpha * p * A * L^0.5))).^(1/(alpha - 1));

% Oferta

b_promedio_optimo_por_r_movil = NaN(N_r_grid, 1);


for r_movil=linspace(1, 1.2, N_r_grid);
    
    % Función de Montecarlo
    z2 = montecarlo_ingresos(0,sqrt(1/2),N_sim,rho,z1);
    
    for i=1:N_sim
        b_grid = linspace(-z2(i)/r_movil, z1, N_b_grid+2)';
        b_grid = b_grid(2:end-1);
        max_b = find(...
            log(z1 - b_grid) + beta * log(b_grid*r_movil + z2(i)) == max(log(z1 - b_grid) + beta * log(b_grid*r_movil + z2(i)))...
            );
        if length(max_b) ~= 1  % Si hay más de un máximo
            b_optimo_por_sim(i) = b_grid(max_b(1));  % toma el primero (supuesto)
        else
            b_optimo_por_sim(i) = b_grid(max_b);
        end
    end
    
    b_promedio_optimo_por_r_movil(sum(isnan(b_promedio_optimo_por_r_movil)==0)+1) = mean(b_optimo_por_sim);

end


% Plot
% Demanda
plot(r_grid, Kt, 'Color', 'black');
set(gca, 'Box', 'off');
hold on; 

ylab = ylabel("$$r$$", 'FontSize', 15, 'Position',[1 4.2515], 'Interpreter', 'latex');
set(ylab,'rotation',0,'VerticalAlignment','bottom')

% Oferta
plot(r_grid, b_promedio_optimo_por_r_movil, 'Color', 'black');

% Cambiar posición
xlabel("$$K_t, b_t$$", 'FontSize', 15, 'Position', [1.21 2.9], 'Interpreter', 'latex') 

hold off;



%% Pregunta 3.1

clear all; close all;
rng(73);  % Seed

P = [...
    0.5	0.5 0 0 0 0;...
    0   0   1 0 0 0;...
    1/3 0   0 1/3 1/3 0;...
    0   0   0 1/2 1/2 0;...
    0   0.1 0 0   0   0.9;...
    0   0   0 0.8 0.2 0;...
    ]

periodos_primera_cadena = 200;
periodos_segunda_cadena = 600;
periodos_tercera_cadena = 1500;


estado_principio = 1;
cadena_1 = Markov(P, estado_principio, periodos_primera_cadena);
cadena_2 = Markov(P, estado_principio, periodos_segunda_cadena);
cadena_3 = Markov(P, estado_principio, periodos_tercera_cadena);

plot((1:periodos_primera_cadena), cadena_1, "r")
hold on;
ylab = ylabel("$Estado$", 'Interpreter', 'latex');
xlabel("$$Pasos$$", 'Interpreter', 'latex') 
hold off;

plot((1:periodos_segunda_cadena), cadena_2, "r")
hold on;
ylab = ylabel("$Estado$", 'Interpreter', 'latex');
xlabel("$$Pasos$$", 'Interpreter', 'latex') 
hold off;

plot((1:periodos_tercera_cadena), cadena_3, "r")
hold on;
ylab = ylabel("$Estado$", 'Interpreter', 'latex');
xlabel("$$Pasos$$", 'Interpreter', 'latex') 
hold off;


%% Pregunta 3.2

% Análisis de convergencia

tamanho_periodos = 1000;

N_iteraciones = 100;
ultimas_posiciones_100 = Markov_iteration_last(P, estado_principio, tamanho_periodos, N_iteraciones);
ultimas_posiciones_100 
ultimas_posiciones_100_normalizadas = normalizadora(...
    unique(ultimas_posiciones_100), ultimas_posiciones_100...
    )

N_iteraciones = 500;
ultimas_posiciones_500 = Markov_iteration_last(P, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_500), ultimas_posiciones_500...
    )

N_iteraciones = 1000;
ultimas_posiciones_1000 = Markov_iteration_last(P, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_1000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_1000), ultimas_posiciones_1000...
    )

N_iteraciones = 2500;
ultimas_posiciones_2500 = Markov_iteration_last(P, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_2500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_2500), ultimas_posiciones_2500...
    )

N_iteraciones =5000;
ultimas_posiciones_5000 = Markov_iteration_last(P, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_5000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_5000), ultimas_posiciones_5000...
    )

plot(ultimas_posiciones_100_normalizadas(:,1), ultimas_posiciones_100_normalizadas(:,2), "r")
hold on;
plot(ultimas_posiciones_500_normalizadas(:,1), ultimas_posiciones_500_normalizadas(:,2), "g")
plot(ultimas_posiciones_1000_normalizadas(:,1), ultimas_posiciones_1000_normalizadas(:,2), "b")
plot(ultimas_posiciones_2500_normalizadas(:,1), ultimas_posiciones_2500_normalizadas(:,2), "c")
plot(ultimas_posiciones_5000_normalizadas(:,1), ultimas_posiciones_5000_normalizadas(:,2), "m")
legend(...
    '100 iteraciones',...
    '500 iteraciones',...
    '1000 iteraciones',...
    '2500 iteraciones',...
    '5000 iteraciones',...
    'Location','southeast' ...
    )
ylab = ylabel("$Frecuencia$", 'Interpreter', 'latex');
xlabel("$$Estado$$", 'Interpreter', 'latex') 

hold off;

% histogram(ultimas_posiciones_100,'Normalization','probability')
% hold on;
% histogram(ultimas_posiciones_500,'Normalization','probability')
% histogram( ultimas_posiciones_1000,'Normalization','probability')
% histogram( ultimas_posiciones_2500,'Normalization','probability')
% histogram( ultimas_posiciones_5000,'Normalization','probability')
% legend(...
%     '100 iteraciones',...
%     '500 iteraciones',...
%     '1000 iteraciones',...
%     '2500 iteraciones',...
%     '5000 iteraciones'...
%     )
% hold off;



%% Pregunta 3.3


% elevar_p(P, 2) == P*P
% --> true

% elevar_p(P, 2)*P == P*P*P
% --> true
% elevar_p(P, 3) == P*P*P
% --> false ==> las diferencias se deben a aproximaciones intermedias

% elevar_p(P, 2)*P*P == P*P*P*P
% --> true
% elevar_p(P, 4) == P*P*P*P
% --> false ==> las diferencias se deben a aproximaciones intermedias

P_elevado_10 = elevar_p(P, 10)
P_elevado_50 = elevar_p(P, 50)
P_elevado_100 = elevar_p(P, 100)
P_elevado_1000 = elevar_p(P, 1000)


%% Pregunta 3.4


figure;

subplot(2,2,1);
hold on;

% P^10

N_iteraciones = 100;
ultimas_posiciones_100 = Markov_iteration_last(P_elevado_10, estado_principio, tamanho_periodos, N_iteraciones);
ultimas_posiciones_100 
ultimas_posiciones_100_normalizadas = normalizadora(...
    unique(ultimas_posiciones_100), ultimas_posiciones_100...
    )

N_iteraciones = 500;
ultimas_posiciones_500 = Markov_iteration_last(P_elevado_10, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_500), ultimas_posiciones_500...
    )

N_iteraciones = 1000;
ultimas_posiciones_1000 = Markov_iteration_last(P_elevado_10, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_1000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_1000), ultimas_posiciones_1000...
    )

N_iteraciones = 2500;
ultimas_posiciones_2500 = Markov_iteration_last(P_elevado_10, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_2500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_2500), ultimas_posiciones_2500...
    )

N_iteraciones =5000;
ultimas_posiciones_5000 = Markov_iteration_last(P_elevado_10, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_5000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_5000), ultimas_posiciones_5000...
    )

plot(ultimas_posiciones_100_normalizadas(:,1), ultimas_posiciones_100_normalizadas(:,2), "r")
hold on;
plot(ultimas_posiciones_500_normalizadas(:,1), ultimas_posiciones_500_normalizadas(:,2), "g")
plot(ultimas_posiciones_1000_normalizadas(:,1), ultimas_posiciones_1000_normalizadas(:,2), "b")
plot(ultimas_posiciones_2500_normalizadas(:,1), ultimas_posiciones_2500_normalizadas(:,2), "c")
plot(ultimas_posiciones_5000_normalizadas(:,1), ultimas_posiciones_5000_normalizadas(:,2), "m")
legend(...
    '100 iteraciones',...
    '500 iteraciones',...
    '1000 iteraciones',...
    '2500 iteraciones',...
    '5000 iteraciones',...
    'Location','southeast' ...
    )
ylab = ylabel("$Estado$", 'Interpreter', 'latex');
xlabel("$$Pasos$$", 'Interpreter', 'latex') 
title("$P^{10}$",  'Interpreter', 'latex')
hold off;

subplot(2,2,2);
hold on;

% P^50

N_iteraciones = 100;
ultimas_posiciones_100 = Markov_iteration_last(P_elevado_50, estado_principio, tamanho_periodos, N_iteraciones);
ultimas_posiciones_100 
ultimas_posiciones_100_normalizadas = normalizadora(...
    unique(ultimas_posiciones_100), ultimas_posiciones_100...
    )

N_iteraciones = 500;
ultimas_posiciones_500 = Markov_iteration_last(P_elevado_50, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_500), ultimas_posiciones_500...
    )

N_iteraciones = 1000;
ultimas_posiciones_1000 = Markov_iteration_last(P_elevado_50, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_1000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_1000), ultimas_posiciones_1000...
    )

N_iteraciones = 2500;
ultimas_posiciones_2500 = Markov_iteration_last(P_elevado_50, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_2500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_2500), ultimas_posiciones_2500...
    )

N_iteraciones =5000;
ultimas_posiciones_5000 = Markov_iteration_last(P_elevado_50, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_5000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_5000), ultimas_posiciones_5000...
    )

plot(ultimas_posiciones_100_normalizadas(:,1), ultimas_posiciones_100_normalizadas(:,2), "r")
hold on;
plot(ultimas_posiciones_500_normalizadas(:,1), ultimas_posiciones_500_normalizadas(:,2), "g")
plot(ultimas_posiciones_1000_normalizadas(:,1), ultimas_posiciones_1000_normalizadas(:,2), "b")
plot(ultimas_posiciones_2500_normalizadas(:,1), ultimas_posiciones_2500_normalizadas(:,2), "c")
plot(ultimas_posiciones_5000_normalizadas(:,1), ultimas_posiciones_5000_normalizadas(:,2), "m")
legend(...
    '100 iteraciones',...
    '500 iteraciones',...
    '1000 iteraciones',...
    '2500 iteraciones',...
    '5000 iteraciones',...
    'Location','southeast' ...
    )
ylab = ylabel("$Estado$", 'Interpreter', 'latex');
xlabel("$$Pasos$$", 'Interpreter', 'latex') 
title("$P^{50}$",  'Interpreter', 'latex')
hold off;

subplot(2,2,3);
hold on;

% P^100

N_iteraciones = 100;
ultimas_posiciones_100 = Markov_iteration_last(P_elevado_100, estado_principio, tamanho_periodos, N_iteraciones);
ultimas_posiciones_100 
ultimas_posiciones_100_normalizadas = normalizadora(...
    unique(ultimas_posiciones_100), ultimas_posiciones_100...
    )

N_iteraciones = 500;
ultimas_posiciones_500 = Markov_iteration_last(P_elevado_100, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_500), ultimas_posiciones_500...
    )

N_iteraciones = 1000;
ultimas_posiciones_1000 = Markov_iteration_last(P_elevado_100, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_1000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_1000), ultimas_posiciones_1000...
    )

N_iteraciones = 2500;
ultimas_posiciones_2500 = Markov_iteration_last(P_elevado_100, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_2500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_2500), ultimas_posiciones_2500...
    )

N_iteraciones =5000;
ultimas_posiciones_5000 = Markov_iteration_last(P_elevado_100, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_5000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_5000), ultimas_posiciones_5000...
    )

plot(ultimas_posiciones_100_normalizadas(:,1), ultimas_posiciones_100_normalizadas(:,2), "r")
hold on;
plot(ultimas_posiciones_500_normalizadas(:,1), ultimas_posiciones_500_normalizadas(:,2), "g")
plot(ultimas_posiciones_1000_normalizadas(:,1), ultimas_posiciones_1000_normalizadas(:,2), "b")
plot(ultimas_posiciones_2500_normalizadas(:,1), ultimas_posiciones_2500_normalizadas(:,2), "c")
plot(ultimas_posiciones_5000_normalizadas(:,1), ultimas_posiciones_5000_normalizadas(:,2), "m")
legend(...
    '100 iteraciones',...
    '500 iteraciones',...
    '1000 iteraciones',...
    '2500 iteraciones',...
    '5000 iteraciones',...
    'Location','southeast' ...
    )
ylab = ylabel("$Estado$", 'Interpreter', 'latex');
xlabel("$$Pasos$$", 'Interpreter', 'latex') 
title("$P^{100}$",  'Interpreter', 'latex')
hold off;


subplot(2,2,4);
hold on;

% P^1000



N_iteraciones = 100;
ultimas_posiciones_100 = Markov_iteration_last(P_elevado_1000, estado_principio, tamanho_periodos, N_iteraciones);
ultimas_posiciones_100 
ultimas_posiciones_100_normalizadas = normalizadora(...
    unique(ultimas_posiciones_100), ultimas_posiciones_100...
    )

N_iteraciones = 500;
ultimas_posiciones_500 = Markov_iteration_last(P_elevado_1000, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_500), ultimas_posiciones_500...
    )

N_iteraciones = 1000;
ultimas_posiciones_1000 = Markov_iteration_last(P_elevado_1000, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_1000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_1000), ultimas_posiciones_1000...
    )

N_iteraciones = 2500;
ultimas_posiciones_2500 = Markov_iteration_last(P_elevado_1000, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_2500_normalizadas = normalizadora(...
    unique(ultimas_posiciones_2500), ultimas_posiciones_2500...
    )

N_iteraciones =5000;
ultimas_posiciones_5000 = Markov_iteration_last(P_elevado_1000, estado_principio, tamanho_periodos, N_iteraciones)
ultimas_posiciones_5000_normalizadas = normalizadora(...
    unique(ultimas_posiciones_5000), ultimas_posiciones_5000...
    )

plot(ultimas_posiciones_100_normalizadas(:,1), ultimas_posiciones_100_normalizadas(:,2), "r")
hold on;
plot(ultimas_posiciones_500_normalizadas(:,1), ultimas_posiciones_500_normalizadas(:,2), "g")
plot(ultimas_posiciones_1000_normalizadas(:,1), ultimas_posiciones_1000_normalizadas(:,2), "b")
plot(ultimas_posiciones_2500_normalizadas(:,1), ultimas_posiciones_2500_normalizadas(:,2), "c")
plot(ultimas_posiciones_5000_normalizadas(:,1), ultimas_posiciones_5000_normalizadas(:,2), "m")
legend(...
    '100 iteraciones',...
    '500 iteraciones',...
    '1000 iteraciones',...
    '2500 iteraciones',...
    '5000 iteraciones',...
    'Location','southeast' ...
    )
ylab = ylabel("$Estado$", 'Interpreter', 'latex');
xlabel("$$Pasos$$", 'Interpreter', 'latex') 
title("$P^{1000}$",  'Interpreter', 'latex')
hold off;

