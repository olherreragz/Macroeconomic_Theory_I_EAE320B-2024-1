function datos_normalizados = normalizadora(estados, vector)

    conteo = tabulate(categorical(vector))
    datos_normalizados = zeros(length(estados), 2);
    datos_normalizados(:,1) = estados;
    
    for i = 1:length(estados)
        datos_normalizados(i,2) = conteo{i,3};
    end    
    
    datos_normalizados(:,2) = datos_normalizados(:,2)/100;
    
end
