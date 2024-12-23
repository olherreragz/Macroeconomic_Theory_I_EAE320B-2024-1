function vector_de_ultimos_estados = Markov_iteration_last(Matriz, estado_actual, T, iterations)

    vector_de_ultimos_estados = NaN(iterations, 1);
    
    for i=1:iterations
        vector_de_ultimos_estados(i) = Markov_last(Matriz, estado_actual, T);
    end
    
end
