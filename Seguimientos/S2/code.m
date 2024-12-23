%% Pregunta 1

rng(73);  % Seed

% Generación de la serie

T = 53;
t = 1:T;
t = t';
estado_laboral = zeros(T,1);

estado = 2;
estado_laboral(1,1) = estado;


for i = 1:T-1  % iter T-1 veces
    if estado == 1
        simulacion = rand;
        if simulacion<=0.87
            estado = 1;
        else
            estado = 2;
        end
    else 
        simulacion = rand;
        if simulacion<=0.6
            estado = 1;
        else
            estado = 2;
        end
    end
    estado_laboral(i+1,1) = estado;
end


% Desplegar los 10 primeros estados y los 10 últimos estados (con ID)
display_serie = [t estado_laboral];
[display_serie(1:10,:); display_serie(length(display_serie)-9:end,:)]


%% Pregunta 2

% Gráfico de trayectoria

figure;
plot(display_serie(:,1),display_serie(:,2),'-gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]...
    )
hold on;
xlim([0 54]);
yticks([1 1.5 2]);
title('Trayectoria para el desempleo de la economía.', 'FontSize', 22, 'Position',[54/2 2.05])
ylabel('Estado de desempleo', 'Position',[-3 1.5],'Fontsize',20)
xlabel('Periodo','Position',[26 0.9],'Fontsize',20)
xticklabel = get(gca,'XTickLabel');
set(gca,'XTickLabel',xticklabel,'fontsize',18)
xticks([0 5 10 15 20 25 30 35 40 45 50]);
hold off;
% Se corren 60.5 pixeles a la izquierda a mano en el código del svg


%% Pregunta 3

% Cálculo de porcentaje
sum(estado_laboral > 1)/length(estado_laboral)


