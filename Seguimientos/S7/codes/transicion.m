function [new_state] = transicion(P, estado_actual)
    numero_aleatorio = rand;
    suma_acumulada = 0;
    for i = 1:size(P,2)
        if P(estado_actual, i) == 0
            continue;
        else 
            suma_acumulada = suma_acumulada + P(estado_actual,i);
            if numero_aleatorio <= suma_acumulada
                new_state = i;
                break
            end
        end
    end
end