function valor_a_devolver = Markov_last(Matriz, estado_actual, T)

    cadena_simulada = Markov(Matriz, estado_actual, T);
    valor_a_devolver = cadena_simulada(end);
    
end
