function V_norm = normalize_vectors_safe(V)

    norms = vecnorm(V, 2, 2);
    norms(norms == 0) = 1;

    V_norm = V ./ norms;

end