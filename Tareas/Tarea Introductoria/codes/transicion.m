function nuevo_estado = transicion(Matriz, estado_actual)

    dimensiones_matriz = size(Matriz);
    nro_estados_posibles = dimensiones_matriz(2);

    rango_de_simulacion_asociado_al_estado = NaN(nro_estados_posibles,1);
%     simulacion = rand
    simulacion = rand;
    
   
    for i=1:nro_estados_posibles
%         if Matriz(estado_actual, i) == 1.01;
%             rango_de_simulacion_asociado_al_estado(i) = 0;
%         elseif i>1 & sum(Matriz(estado_actual, 1:i-1))==1  % Probablemente tiré error de índice para i=1
%             rango_de_simulacion_asociado_al_estado(i) = 1.01;
        if i>1 & sum(Matriz(estado_actual, 1:i-1))==1  % Probablemente tiré error de índice para i=1
            rango_de_simulacion_asociado_al_estado(i) = 1.01;
        else
            rango_de_simulacion_asociado_al_estado(i) = sum(Matriz(estado_actual, 1:i));
        end
    end
    
%     [(1:nro_estados_posibles)' rango_de_simulacion_asociado_al_estado]
    
%     (rango_de_simulacion_asociado_al_estado <= simulacion)
%     nuevo_estado = sum(rango_de_simulacion_asociado_al_estado <= simulacion) + 1
    nuevo_estado = sum(rango_de_simulacion_asociado_al_estado <= simulacion) + 1;

    
end
