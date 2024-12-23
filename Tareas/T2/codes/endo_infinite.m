function [I_policy] = endo_infinite (gridK, beta, alpha, a, A, Delta, gamma, tol, Kmin, Kmax)

    % Condición inicial KT+1 = 0 + mov. de capital:
    IT = -gridK*(1-Delta)
    
    It = (1/(2*gamma)) * (beta*(...
        alpha*a*A*gridK.^(alpha - 1) -...
        2*alpha*A^2*gridK.^(2*alpha - 1) + ...
        (1-Delta)*(1+2*gamma*IT)...
        )...
        - 1);    
    % No censuramos a valores positivos porque podría "desinvertir"

    
    % Despejando Kt ∀t:
    % Kt = (Kt+1 - It)/1(1-Delta)
    KT = (gridK - IT)/(1-Delta);

    % Interpolación para despejar Policy en función de KT's de nuestro
    % interés
    It = interp1(KT,IT,gridK,[],'extrap'); 

    % CHECK que KT+1 amarrado esté en el rango permitido
    % Kt+1 = (1-Delta)*Kt + It
    K_T_plus_1_amarrado = (1-Delta)*gridK + It;

    K_T_plus_1_amarrado(K_T_plus_1_amarrado<Kmin) = Kmin
    
    % Una vez que forzamos a que no sea menor al mínimo, calculamos
    % la inversión que nos asegura no pasarnos en T+1
    It = K_T_plus_1_amarrado - (1-Delta)*gridK

    I_policy(:,2) = It;

    % Ahora mi It pasa a ser mi It+1,
    % y mi Kt asociado a It, pasa a ser el Kt+1 asociado a It+1:
    It = (1/(2*gamma)) * (beta*(...
        alpha*a*A*gridK.^(alpha - 1) -...
        2*alpha*A^2*gridK.^(2*alpha - 1) + ...
        (1-Delta)*(1+2*gamma*It)...
        )...
        - 1);
    
    % Despejando Kt:
    KT = (gridK - It)/(1-Delta);  % Acá se reemplaza It por IT

    % Interpolación para despejar Policy
    It = interp1(KT,IT,gridK,[],'extrap');  % Acá se reemplaza It por IT
    
    % CHECK Kt+1
    K_T_plus_1_amarrado = (1-Delta)*gridK + It;


    K_T_plus_1_amarrado(K_T_plus_1_amarrado<Kmin) = Kmin
    
    It = K_T_plus_1_amarrado - (1-Delta)*gridK

    I_policy(:,3) = It;
    
    % Iterando:
    D = 100000;
    i = 4;
    while D > tol
        disp(i)
        It = (1/(2*gamma)) * (beta*(...
            alpha*a*A*gridK.^(alpha - 1) -...
            2*alpha*A^2*gridK.^(2*alpha - 1) + ...
            (1-Delta)*(1+2*gamma*It)...
            )...
            - 1);
        KT = (gridK - It)/(1-Delta);
        It = interp1(KT,It,gridK,[],'extrap');  % Acá se reemplaza It por IT
        
        % CHECK Kt+1
        K_T_plus_1_amarrado = (1-Delta)*gridK + It;

        K_T_plus_1_amarrado(K_T_plus_1_amarrado<Kmin) = Kmin

        It = K_T_plus_1_amarrado - (1-Delta)*gridK

        I_policy(:,i) = It;
        
        D = max(abs(I_policy(:,i) - I_policy(:,i-1)));
        % max(abs(CC(:,:,i) - CC(:,:,i-1))) toma la máxima distancia entre 
        % las policies para los dos posibles estados

        % max(abs(I_policy(:,i) - I_policy(:,i-1))) nos asegura parar cuando
        % todas las policies (para todos los capitales) hayan convergido

        i = i+1;
        disp(D)
    end    

end