function [CC] = endo_infinite (gridA, N_A, gridZ, N_Z, w, P, beta, sigma, R, phi)
    theta = w.*gridZ'; % Senda de ingresos para el estado alto y bajo 
    % Condición inicial A = C
    % <==> Aj,T+1 = 0  ==>  C_j = A*R + w*Z_j
    
%     repmat(w.*gridZ',1,size(gridA, 2)) + repmat(gridA*R,2,1) <-- ∀A,
%            genero posible CT según realización de Z
%            (Riqueza líquida en T)

    CT = repmat(theta,1,N_A) + repmat(gridA,N_Z,1).*R;
%   CT == repmat(w.*gridZ',1,size(gridA, 2)) + repmat(gridA*R,2,1)
    CT = CT.*(CT > 0);
    % Tenemos que las filas son los estados, las columnas son las grillas (Ahorros)
    % Para agregar la iteración agregamos una tercera dimención y nos queda
    % que los consumos optimos CC serian los siguientes
    CC(:,:,1) = CT;  % Al principio, al generar CC(...), CC es una matriz
    % CC(estado_naturaleza, grilla_de_estados, iteracion)
    
    E = P*(CT.^(-sigma)); 
    C = ((beta*R).^(-inv(sigma)))*(E.^(-inv(sigma)));
    C = C.*(C > 0);
    % M = AT-1 + w*Zj dado AT
    M = repmat(gridA,N_Z,1) + C;
    A = (M - repmat(theta,1,N_A))./R;
    % A = AT-1  (stock de activos corrientes)
    for j=1:N_Z
    % nos da el consumo iterpolado en base a la grilla original.
    % no necesariamente nuestro A de antes es igual al de la grilla.
        CT(j,:) = interp1(A(j,:),C(j,:),gridA,[],'extrap'); 
    end
    % pisamos CT para usarlo como auxiliar en cada iteración
    
    % Dado un Ct, podemos despejar con qué ahorro inicial At es consistente
    % ese consumo, dado un At+1
    AT = repmat(theta,1,N_A) + repmat(gridA,N_Z,1).*R - CT;
    % Forzamos que ese "AT consistente" cumpla la restricción
    AT(AT<-phi) = -phi; 
    
    % Ajustamos el consumo si es que la restricción está activa
    % para algún consumo CT encontrado
    CT = repmat(theta,1,N_A) + repmat(gridA,N_Z,1).*R - AT;
    CC(:,:,2) = CT;
    
    % Iterando:
    D = 100000; i = 3;
    while D >.0001
        disp(i)
    E = P*(CT.^(-sigma))
    C = ((beta*R).^(-inv(sigma)))*(E.^(-inv(sigma)));
    C = C.*(C > 0);
    M = repmat(gridA,N_Z,1) + C;
    A = (M - repmat(theta,1,N_A))./R;
    for j=1:N_Z
        CT(j,:) = interp1(A(j,:),C(j,:),gridA,[],'extrap'); 
    end
    
    AT = repmat(theta,1,N_A) + repmat(gridA,N_Z,1).*R - CT;
    AT(AT<-phi) = -phi; 
    CT = repmat(theta,1,N_A) + repmat(gridA,N_Z,1).*R - AT;
    
    CC(:,:,i) = CT;
    D = max(max(abs(CC(:,:,i) - CC(:,:,i-1))));
    % max(abs(CC(:,:,i) - CC(:,:,i-1))) toma la máxima distancia entre 
    % las policies para los dos posibles estados

    % max(max(abs(CC(:,:,i) - CC(:,:,i-1)))) nos asegura parar cuando
    % todas las policies, para todos los estados y activos, hayan convergido

    i = i+1;
    disp(D)
    end
    
end
