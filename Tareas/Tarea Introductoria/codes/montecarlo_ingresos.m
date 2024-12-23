function z2 = montecarlo_ingresos(mean_z2, sd_z2, N_simulaciones,rho,z1)

    epsilon = normrnd(mean_z2,sd_z2,N_simulaciones,1);
    z2 = rho * ones(length(epsilon), 1) * sqrt(z1) + epsilon;

end
