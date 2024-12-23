function utils= udd_cobb_douglas(alpha,X1,X2)

    % x1^\alpha*x2^(1-\alpha)
    utils = X1.^alpha.*X2.^(1-alpha)

end