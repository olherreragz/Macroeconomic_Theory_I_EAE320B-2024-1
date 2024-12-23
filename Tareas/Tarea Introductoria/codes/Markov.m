function cadena_de_markov = Markov(Matriz, estado_actual, T)

    cadena_de_markov = NaN(T,1);
    cadena_de_markov(1) = estado_actual;
    
    for i=2:T;
        cadena_de_markov(i) = transicion(Matriz,cadena_de_markov(i-1));
    end
    
end
