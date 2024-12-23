function matriz_final = elevar_p(Matriz, n)

%     Caso base
    if n == 0
        matriz_final = eye(size(Matriz));
    elseif n == 1
        matriz_final = Matriz;
%     if n == 1
%         matriz_final = Matriz;
    else
%         matriz_intermedia = elevar_p(Matriz, n-1);
%         matriz_final = Matriz * matriz_intermedia;
        matriz_final = Matriz * elevar_p(Matriz, n-1);;

end


