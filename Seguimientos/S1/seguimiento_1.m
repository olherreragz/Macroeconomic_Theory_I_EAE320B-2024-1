%% Seguimiento 1 Oscar Herrera

clc; close all; clear;
rng(73);  % Seed
n = 418;
w = 4;

% Salas ∼ U(30, 100)
salas = 30 + (100-30)*rand(n,w);
display_salas = [(1:418)' salas];
display_salas(1:10,:)


% Estudio ∼ N(0, 3600)
% raw_estudio = normrnd(0,3600,n,w);
raw_estudio = normrnd(0,60,n,w);
display_raw_estudio = [(1:418)' raw_estudio];
display_raw_estudio(1:10,:)


% Viaje ∼ U(1, 600)
viaje = 1 + (600-1)*rand(n,w);
display_viaje = [(1:418)' viaje];
display_viaje(1:10,:)


%% 1. Crecion de M_DECANO_LOOP

M_DECANO_LOOP= zeros(n,w);
for i= 1:n
    for j = 1:w
        if raw_estudio(i,j)<0
            M_DECANO_LOOP(i,j) = salas(i,j) + 0 + viaje(i,j);
        else
            M_DECANO_LOOP(i,j) = salas(i,j) + raw_estudio(i,j) + viaje(i,j);
        end
    end
end


%% 2. Creacion de M_DECANO_MATRICIAL

estudio = (1-(raw_estudio<0)).*raw_estudio;
M_DECANO_MATRICIAL = salas + estudio + viaje
M_DECANO_MATRICIAL


%% 3. Chequeo de igualdad entre matrices

chequeo = sum(M_DECANO_MATRICIAL == M_DECANO_LOOP)

if chequeo(1,1)==n & chequeo(1,2)==n & chequeo(1,3)==n & chequeo(1,4)==n
    display("Efectivamente, son iguales")
else
    display("no son iguales")
end

display_M_DECANO_LOOP = [(1:418)' M_DECANO_LOOP];
display_M_DECANO_LOOP(1:10,:)


%% 4. Identificando al estudiante que más tiempo usó

% Concatenamos un vector columna de ids para poder identificar el alumno
MATRIZ_SALAS = [(1:418)' salas];
MATRIZ_ESTUDIO = [(1:418)' estudio];
MATRIZ_VIAJE = [(1:418)' viaje];
MATRIZ_TOTAL = [(1:418)' M_DECANO_MATRICIAL];

% Estudiante que más tiempo gastó en cambio de salas
max_tiempo_sala_por_semana = max(salas);
max_tiempo_sala_w1 = max_tiempo_sala_por_semana(1,1);
id_max_tiempo_sala_w1  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,2) == max_tiempo_sala_w1))
display("El estudiante que más tiempo usó en cambio de salas en la semana 1 es el estudiante nro. " + id_max_tiempo_sala_w1)

max_tiempo_sala_w2 = max_tiempo_sala_por_semana(1,2);
id_max_tiempo_sala_w2  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,3) == max_tiempo_sala_w2))
display("El estudiante que más tiempo usó en cambio de salas en la semana 2 es el estudiante nro. " + id_max_tiempo_sala_w2)

max_tiempo_sala_w3 = max_tiempo_sala_por_semana(1,3);
id_max_tiempo_sala_w3  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,4) == max_tiempo_sala_w3))
display("El estudiante que más tiempo usó en cambio de salas en la semana 3 es el estudiante nro. " + id_max_tiempo_sala_w3)

max_tiempo_sala_w4 = max_tiempo_sala_por_semana(1,4);
id_max_tiempo_sala_w4  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,5) == max_tiempo_sala_w4))
display("El estudiante que más tiempo usó en cambio de salas en la semana 4 es el estudiante nro. " + id_max_tiempo_sala_w4)

maximos_sala_por_semana = [...
    [1 2 3 4]'...
    [id_max_tiempo_sala_w1 id_max_tiempo_sala_w2 id_max_tiempo_sala_w3 id_max_tiempo_sala_w4]'...
    [max_tiempo_sala_w1 max_tiempo_sala_w2 max_tiempo_sala_w3 max_tiempo_sala_w4]'...
    ]



% Estudiante que más tiempo gastó en estudio
max_tiempo_estudio_por_semana = max(estudio);
max_tiempo_estudio_w1 = max_tiempo_estudio_por_semana(1,1);
id_max_tiempo_estudio_w1  = max(MATRIZ_ESTUDIO(:,1).*(MATRIZ_ESTUDIO(:,2) == max_tiempo_estudio_w1))
display("El estudiante que más tiempo usó en estudio en la semana 1 es el estudiante nro. " + id_max_tiempo_estudio_w1)

max_tiempo_estudio_w2 = max_tiempo_estudio_por_semana(1,2);
id_max_tiempo_estudio_w2  = max(MATRIZ_ESTUDIO(:,1).*(MATRIZ_ESTUDIO(:,3) == max_tiempo_estudio_w2))
display("El estudiante que más tiempo usó en estudio en la semana 2 es el estudiante nro. " + id_max_tiempo_estudio_w2)

max_tiempo_estudio_w3 = max_tiempo_estudio_por_semana(1,3);
id_max_tiempo_estudio_w3  = max(MATRIZ_ESTUDIO(:,1).*(MATRIZ_ESTUDIO(:,4) == max_tiempo_estudio_w3))
display("El estudiante que más tiempo usó en estudio en la semana 3 es el estudiante nro. " + id_max_tiempo_estudio_w3)

max_tiempo_estudio_w4 = max_tiempo_estudio_por_semana(1,4);
id_max_tiempo_estudio_w4  = max(MATRIZ_ESTUDIO(:,1).*(MATRIZ_ESTUDIO(:,5) == max_tiempo_estudio_w4))
display("El estudiante que más tiempo usó en estudio en la semana 4 es el estudiante nro. " + id_max_tiempo_estudio_w4)

maximos_estudio_por_semana = [...
    [1 2 3 4]'...
    [id_max_tiempo_estudio_w1 id_max_tiempo_estudio_w2 id_max_tiempo_estudio_w3 id_max_tiempo_estudio_w4]'...
    [max_tiempo_estudio_w1 max_tiempo_estudio_w2 max_tiempo_estudio_w3 max_tiempo_estudio_w4]'...
    ]



% Estudiante que más tiempo gastó en viaje
max_tiempo_viaje_por_semana = max(viaje);
max_tiempo_viaje_w1 = max_tiempo_viaje_por_semana(1,1);
id_max_tiempo_viaje_w1  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,2) == max_tiempo_viaje_w1))
display("El estudiante que más tiempo usó en viaje en la semana 1 es el estudiante nro. " + id_max_tiempo_viaje_w1)

max_tiempo_viaje_w2 = max_tiempo_viaje_por_semana(1,2);
id_max_tiempo_viaje_w2  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,3) == max_tiempo_viaje_w2))
display("El estudiante que más tiempo usó en viaje en la semana 2 es el estudiante nro. " + id_max_tiempo_viaje_w2)

max_tiempo_viaje_w3 = max_tiempo_viaje_por_semana(1,3);
id_max_tiempo_viaje_w3  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,4) == max_tiempo_viaje_w3))
display("El estudiante que más tiempo usó en viaje en la semana 3 es el estudiante nro. " + id_max_tiempo_viaje_w3)

max_tiempo_viaje_w4 = max_tiempo_viaje_por_semana(1,4);
id_max_tiempo_viaje_w4  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,5) == max_tiempo_viaje_w4))
display("El estudiante que más tiempo usó en viaje en la semana 4 es el estudiante nro. " + id_max_tiempo_viaje_w4)

maximos_viaje_por_semana = [...
    [1 2 3 4]'...
    [id_max_tiempo_viaje_w1 id_max_tiempo_viaje_w2 id_max_tiempo_viaje_w3 id_max_tiempo_viaje_w4]'...
    [max_tiempo_viaje_w1 max_tiempo_viaje_w2 max_tiempo_viaje_w3 max_tiempo_viaje_w4]'...
    ]



% Estudiante que más tiempo gastó en total
max_tiempo_total_por_semana = max(M_DECANO_MATRICIAL);
max_tiempo_total_w1 = max_tiempo_total_por_semana(1,1);
id_max_tiempo_total_w1  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,2) == max_tiempo_total_w1))
display("El estudiante que más tiempo usó en total en la semana 1 es el estudiante nro. " + id_max_tiempo_total_w1)

max_tiempo_total_w2 = max_tiempo_total_por_semana(1,2);
id_max_tiempo_total_w2  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,3) == max_tiempo_total_w2))
display("El estudiante que más tiempo usó en total en la semana 2 es el estudiante nro. " + id_max_tiempo_total_w2)

max_tiempo_total_w3 = max_tiempo_total_por_semana(1,3);
id_max_tiempo_total_w3  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,4) == max_tiempo_total_w3))
display("El estudiante que más tiempo usó en total en la semana 3 es el estudiante nro. " + id_max_tiempo_total_w3)

max_tiempo_total_w4 = max_tiempo_total_por_semana(1,4);
id_max_tiempo_total_w4  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,5) == max_tiempo_total_w4))
display("El estudiante que más tiempo usó en total en la semana 4 es el estudiante nro. " + id_max_tiempo_total_w4)

maximos_total_por_semana = [...
    [1 2 3 4]'...
    [id_max_tiempo_total_w1 id_max_tiempo_total_w2 id_max_tiempo_total_w3 id_max_tiempo_total_w4]'...
    [max_tiempo_total_w1 max_tiempo_total_w2 max_tiempo_total_w3 max_tiempo_total_w4]'...
    ]



%% Suma de todas las semanas por cada variable


MATRIZ_SALAS = [(1:418)' salas];
MATRIZ_ESTUDIO = [(1:418)' estudio];
MATRIZ_VIAJE = [(1:418)' viaje];
MATRIZ_TOTAL = [(1:418)' M_DECANO_MATRICIAL];


MATRIZ_SALAS_MENSUAL = [(1:418)' sum(salas, 2)];
MATRIZ_ESTUDIO_MENSUAL = [(1:418)' sum(estudio, 2)];
MATRIZ_VIAJE_MENSUAL = [(1:418)' sum(viaje, 2)];
MATRIZ_TOTAL_MENSUAL = [(1:418)' sum(M_DECANO_MATRICIAL, 2)];
% format SHORTG

% Estudiante que más tiempo gastó en cambio de salas
max_tiempo_sala_al_mes = max(sum(salas, 2));
id_max_tiempo_sala_mensual = max(MATRIZ_SALAS_MENSUAL(:,1).*(MATRIZ_SALAS_MENSUAL(:,2) == max_tiempo_sala_al_mes))

max_tiempo_estudio_al_mes = max(sum(estudio, 2));
id_max_tiempo_estudio_mensual = max(MATRIZ_ESTUDIO_MENSUAL(:,1).*(MATRIZ_ESTUDIO_MENSUAL(:,2) == max_tiempo_estudio_al_mes))

max_tiempo_viaje_al_mes = max(sum(viaje, 2));
id_max_tiempo_viaje_mensual = max(MATRIZ_VIAJE_MENSUAL(:,1).*(MATRIZ_VIAJE_MENSUAL(:,2) == max_tiempo_viaje_al_mes))

max_tiempo_total_al_mes = max(sum(M_DECANO_MATRICIAL, 2));
id_max_tiempo_total_mensual = max(MATRIZ_TOTAL_MENSUAL(:,1).*(MATRIZ_TOTAL_MENSUAL(:,2) == max_tiempo_total_al_mes))


[...
    ["Tiempo de cambio de salas" "Tiempo de estudio" "Tiempo de viajes" "Tiempo total"]'...
    [id_max_tiempo_sala_mensual id_max_tiempo_estudio_mensual id_max_tiempo_viaje_mensual id_max_tiempo_total_mensual]'...
    [max_tiempo_sala_al_mes max_tiempo_estudio_al_mes max_tiempo_viaje_al_mes max_tiempo_total_al_mes]'...
    ]


%% 5. Tiempos mínimos semanales

% Estudiante que menos tiempo gastó en cambio de salas
min_tiempo_sala_por_semana = min(salas);
min_tiempo_sala_w1 = min_tiempo_sala_por_semana(1,1);
id_min_tiempo_sala_w1  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,2) == min_tiempo_sala_w1))

min_tiempo_sala_w2 = min_tiempo_sala_por_semana(1,2);
id_min_tiempo_sala_w2  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,3) == min_tiempo_sala_w2))

min_tiempo_sala_w3 = min_tiempo_sala_por_semana(1,3);
id_min_tiempo_sala_w3  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,4) == min_tiempo_sala_w3))

min_tiempo_sala_w4 = min_tiempo_sala_por_semana(1,4);
id_min_tiempo_sala_w4  = max(MATRIZ_SALAS(:,1).*(MATRIZ_SALAS(:,5) == min_tiempo_sala_w4))

minimos_sala_por_semana = [...
    [1 2 3 4]'...
    [id_min_tiempo_sala_w1 id_min_tiempo_sala_w2 id_min_tiempo_sala_w3 id_min_tiempo_sala_w4]'...
    [min_tiempo_sala_w1 min_tiempo_sala_w2 min_tiempo_sala_w3 min_tiempo_sala_w4]'...
    ]



% Estudiante que menos tiempo gastó en estudio
min_tiempo_estudio_por_semana = min(estudio);
min_tiempo_estudio_w1 = min_tiempo_estudio_por_semana(1,1);
ids_min_tiempo_estudio_w1  = find(MATRIZ_ESTUDIO(:,2) == min_tiempo_estudio_w1);

length(ids_min_tiempo_estudio_w1)
[ids_min_tiempo_estudio_w1(1:10) ids_min_tiempo_estudio_w1(length(ids_min_tiempo_estudio_w1)-10:end-1)]

min_tiempo_estudio_w2 = min_tiempo_estudio_por_semana(1,2);
ids_min_tiempo_estudio_w2  = find(MATRIZ_ESTUDIO(:,3) == min_tiempo_estudio_w2);

length(ids_min_tiempo_estudio_w2)
[ids_min_tiempo_estudio_w2(1:10) ids_min_tiempo_estudio_w2(length(ids_min_tiempo_estudio_w2)-10:end-1)]

min_tiempo_estudio_w3 = min_tiempo_estudio_por_semana(1,3);
ids_min_tiempo_estudio_w3  = find(MATRIZ_ESTUDIO(:,4) == min_tiempo_estudio_w3);

length(ids_min_tiempo_estudio_w3)
[ids_min_tiempo_estudio_w3(1:10) ids_min_tiempo_estudio_w3(length(ids_min_tiempo_estudio_w3)-10:end-1)]

min_tiempo_estudio_w4 = min_tiempo_estudio_por_semana(1,4);
ids_min_tiempo_estudio_w4  = find(MATRIZ_ESTUDIO(:,5) == min_tiempo_estudio_w4);

length(ids_min_tiempo_estudio_w4)
[ids_min_tiempo_estudio_w4(1:10) ids_min_tiempo_estudio_w4(length(ids_min_tiempo_estudio_w4)-10:end-1)]



% Estudiante que más tiempo gastó en viaje
min_tiempo_viaje_por_semana = min(viaje);
min_tiempo_viaje_w1 = min_tiempo_viaje_por_semana(1,1);
id_min_tiempo_viaje_w1  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,2) == min_tiempo_viaje_w1))

min_tiempo_viaje_w2 = min_tiempo_viaje_por_semana(1,2);
id_min_tiempo_viaje_w2  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,3) == min_tiempo_viaje_w2))

min_tiempo_viaje_w3 = min_tiempo_viaje_por_semana(1,3);
id_min_tiempo_viaje_w3  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,4) == min_tiempo_viaje_w3))

min_tiempo_viaje_w4 = min_tiempo_viaje_por_semana(1,4);
id_min_tiempo_viaje_w4  = max(MATRIZ_VIAJE(:,1).*(MATRIZ_VIAJE(:,5) == min_tiempo_viaje_w4))

minimos_viaje_por_semana = [...
    [1 2 3 4]'...
    [id_min_tiempo_viaje_w1 id_min_tiempo_viaje_w2 id_min_tiempo_viaje_w3 id_min_tiempo_viaje_w4]'...
    [min_tiempo_viaje_w1 min_tiempo_viaje_w2 min_tiempo_viaje_w3 min_tiempo_viaje_w4]'...
    ]



% Estudiante que menos tiempo gastó en total
min_tiempo_total_por_semana = min(M_DECANO_MATRICIAL);
min_tiempo_total_w1 = min_tiempo_total_por_semana(1,1);
id_min_tiempo_total_w1  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,2) == min_tiempo_total_w1))

min_tiempo_total_w2 = min_tiempo_total_por_semana(1,2);
id_min_tiempo_total_w2  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,3) == min_tiempo_total_w2))

min_tiempo_total_w3 = min_tiempo_total_por_semana(1,3);
id_min_tiempo_total_w3  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,4) == min_tiempo_total_w3))

min_tiempo_total_w4 = min_tiempo_total_por_semana(1,4);
id_min_tiempo_total_w4  = max(MATRIZ_TOTAL(:,1).*(MATRIZ_TOTAL(:,5) == min_tiempo_total_w4))

minimos_total_por_semana = [...
    [1 2 3 4]'...
    [id_min_tiempo_total_w1 id_min_tiempo_total_w2 id_min_tiempo_total_w3 id_min_tiempo_total_w4]'...
    [min_tiempo_total_w1 min_tiempo_total_w2 min_tiempo_total_w3 min_tiempo_total_w4]'...
    ]


%% Gráficos
% Máximos

figure;

subplot(2,2,1);
hold on;
bar(maximos_sala_por_semana(:,3),'FaceColor',"#80B3FF",'EdgeColor','none','LineWidth',1.5)
ylim([98 100.5])
xticks([1 2 3 4])
title('Tiempos máximos de cambios de sala')
xlabel('Semana')
box on;
hold off;

subplot(2,2,2);
hold on;
bar(maximos_estudio_por_semana(:,3),'FaceColor',"#80B3FF",'EdgeColor','none','LineWidth',1.5)
xticks([1 2 3 4])
title('Tiempos máximos de estudio')
xlabel('Semana')
box on;
hold off;

subplot(2,2,3);
hold on;
bar(maximos_viaje_por_semana(:,3),'FaceColor',"#80B3FF",'EdgeColor','none','LineWidth',1.5)
ylim([597.5 600.5])
xticks([1 2 3 4])
title('Tiempos máximos de viaje')
xlabel('Semana')
box on;
hold off;

subplot(2,2,4);
hold on;
bar(maximos_total_por_semana(:,3),'FaceColor',"#80B3FF",'EdgeColor','none','LineWidth',1.5)
ylim([750 780])
xticks([1 2 3 4])
title('Tiempos máximos totales')
xlabel('Semana')
box on;
hold off;


% Mínimos

figure;

subplot(2,2,1);
bar(minimos_sala_por_semana(:,3),'FaceColor',"#A1FEC1",'EdgeColor','none','LineWidth',1.5)
ylim([29.5 30.3])
xticks([1 2 3 4])
title('Tiempos mínimos de cambios de sala')
xlabel('Semana')
box on;
hold off;

subplot(2,2,2);
hold on;
scatter([1 2 3 4],[0 0 0 0],100,"square",'filled',"MarkerFaceColor","#A1FEC1")
xticks([1 2 3 4])
xlim([0 5])
title('Tiempos mínimos de estudio')
xlabel('Semana')
box on;
hold off;

subplot(2,2,3);
hold on;
bar(minimos_viaje_por_semana(:,3),'FaceColor',"#A1FEC1",'EdgeColor','none','LineWidth',1.5)
xticks([1 2 3 4])
title('Tiempos mínimos de viaje')
xlabel('Semana')
box on;
hold off;

subplot(2,2,4);
hold on;
bar(minimos_total_por_semana(:,3),'FaceColor',"#A1FEC1",'EdgeColor','none','LineWidth',1.5)
ylim([0 70])
xticks([1 2 3 4])
title('Tiempos mínimos totales')
xlabel('Semana')
box on;
hold off;

