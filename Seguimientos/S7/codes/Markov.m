function [historial_estados] = Markov(P, estado_actual, T)
    historial_estados = zeros(T,1);
    historial_estados(1) = estado_actual;
    for t = 2:T
        historial_estados(t) = transicion(P, estado_actual);
        estado_actual = historial_estados(t);
    end
end
