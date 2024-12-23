function profits = beneficios(p,A,k_grid,alpha,L,r)

    profits = p*A*k_grid.^alpha * L^0.5 - r*k_grid;

end
