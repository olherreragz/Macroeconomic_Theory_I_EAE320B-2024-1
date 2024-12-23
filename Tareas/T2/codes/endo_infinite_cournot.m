function [I_policy] = endo_infinite_cournot (gridK, beta, alpha, a, A, Delta, gamma, tol, Kmin, Kmax)

    grid_cournot_k_i = repmat(gridK', 1, size(gridK, 2));
    grid_cournot_k_j = repmat(gridK, size(gridK, 2), 1);
    
    % Condición inicial KT+1,i = 0,
    % ∀KT+1,j in linspace,
    % & ∀KT,i in grid (linspace);
    
    % Usando la ecuación para el mov. de capital:
    IT = -gridK*(1-Delta);
    % como nuestro KT+1,i inicial es cero
    % y aparte, la ecuación de movimiento no depende de grillas de j
    % IT sólo depende de la variable de estado "personal" en T
  

    % Lo ideal sería guardar sólo una fila 1-by-1000 en I_policy(:,:,1)
    % y guardar matrices en las posiciones de las iteraciones
    % I_policy(:,:,t).
    % Sería algo como:
    % I_policy(:,:,1) = IT;
    % ...
    % I_policy(:,:,t) = It; donde It sería 1000-by-1000 ∀t.
    % Pero no se puede por dimensionalidad (todas tienen que ser 1000x1000).
    % Por lo tanto, se guardará la misma fila repetida 1000 veces en la
    % iteración 1:
    I_policy(:,:,1) = repmat(IT,size(gridK,2), 1);
    IT = repmat(IT,size(gridK,2), 1);
    
    It = (1/(2*gamma)) * (beta*(...
        alpha*a*A*grid_cournot_k_i.^(alpha - 1) -...
        2*alpha*A^2*grid_cournot_k_i.^(2*alpha - 1) - ...
        alpha*A^2*((grid_cournot_k_j.^alpha).*(grid_cournot_k_i.^(alpha - 1))) + ...
        (1-Delta)*(1+2*gamma* IT )...
        )...
        - 1);
    % No censuramos a valores positivos porque podría "desinvertir"

    
    % Despejando Kt ∀t:
    % Kt = (Kt+1 - It)/1(1-Delta)
    
    % Despejando con la ecuación de movimiento el Kt asociado a cada celda:

    KT = (grid_cournot_k_i - It)*(1/(1-Delta));


    % Interpolación para despejar Policy en función de KT's de nuestro
    % interés.
    
    % Fijamos un kj (fijamos columnas de grid_cournot_k_j)
    % then;
    % Ese K_t+1_j tiene asociados distintos K_t_i (KT)
    % que no necesariamente son los mismos valores numéricos que nos interesan
    % para la variable de estado, i.e., no necesariamente son nuestra grilla de interés
    % then;
    % con Kj fijo interpolamos para calcular la respuesta para la grilla
    for j=1:size(grid_cournot_k_j,2)
        It(:,j) = interp1(KT(:,j),It(:,j),grid_cournot_k_i(:,j),[],'extrap'); 
    end

    % CHECK que KT+1 amarrado esté en el rango permitido
    % Kt+1 = (1-Delta)*Kt + It
    K_T_plus_1_amarrado = (1-Delta)*grid_cournot_k_i + It;

    K_T_plus_1_amarrado(K_T_plus_1_amarrado<Kmin) = Kmin;
    
    % Una vez que forzamos a que no sea menor al mínimo, calculamos
    % la inversión que nos asegura no pasarnos en T+1
    It = K_T_plus_1_amarrado - (1-Delta)*grid_cournot_k_i;

    I_policy(:,:,2) = It;
    
    % Ahora mi It pasa a ser mi It+1,
    % y mi Kt asociado a It, pasa a ser el Kt+1 asociado a It+1.
    % Re-iteramos...
    
    It = (1/(2*gamma)) * (beta*(...
        alpha*a*A*grid_cournot_k_i.^(alpha - 1) -...
        2*alpha*A^2*grid_cournot_k_i.^(2*alpha - 1) - ...
        alpha*A^2*((grid_cournot_k_j.^alpha).*(grid_cournot_k_i.^(alpha - 1))) + ...
        (1-Delta)*(1+2*gamma* It )...
        )...
        - 1);

    % Despejando Kt:
    KT = (grid_cournot_k_i - It)*(1/(1-Delta));  % Acá se reemplaza It por IT

    % Interpolación para despejar Policy
    for j=1:size(grid_cournot_k_j,2)
        It(:,j) = interp1(KT(:,j),It(:,j),grid_cournot_k_i(:,j),[],'extrap');   % Acá se reemplaza It por IT
    end
    
    % CHECK Kt+1
    K_T_plus_1_amarrado = (1-Delta)*grid_cournot_k_i + It;

    K_T_plus_1_amarrado(K_T_plus_1_amarrado<Kmin) = Kmin;
    
    It = K_T_plus_1_amarrado - (1-Delta)*grid_cournot_k_i;

    I_policy(:,:,3) = It;
    
    % Iterando:
    D = 100000;
    i = 4;
    max_iter = 10;
%     tol = 0.005;
    while D > tol
        disp(i)
        It = (1/(2*gamma)) * (beta*(...
            alpha*a*A*grid_cournot_k_i.^(alpha - 1) -...
            2*alpha*A^2*grid_cournot_k_i.^(2*alpha - 1) - ...
            alpha*A^2*((grid_cournot_k_j.^alpha).*(grid_cournot_k_i.^(alpha - 1))) + ...
            (1-Delta)*(1+2*gamma* It )...
            )...
            - 1);
        KT = (grid_cournot_k_i - It)*(1/(1-Delta));
        
        for j=1:size(grid_cournot_k_j,2)
            It(:,j) = interp1(KT(:,j),It(:,j),grid_cournot_k_i(:,j),[],'extrap');   % Acá se reemplaza It por IT
        end

        % CHECK Kt+1
        K_T_plus_1_amarrado = (1-Delta)*grid_cournot_k_i + It;

        K_T_plus_1_amarrado(K_T_plus_1_amarrado<Kmin) = Kmin;

        It = K_T_plus_1_amarrado - (1-Delta)*grid_cournot_k_i;

        I_policy(:,:,i) = It;
        
        D = max(max(abs(I_policy(:,:,i) - I_policy(:,:,i-1))));
        % max(abs(CC(:,:,i) - CC(:,:,i-1))) toma la máxima distancia entre 
        % las policies para los dos posibles estados

        % max(abs(I_policy(:,i) - I_policy(:,i-1))) nos asegura parar cuando
        % todas las policies (para todos los capitales) hayan convergido

        i = i+1;
        disp(D)
    end    

end
