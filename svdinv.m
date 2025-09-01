function Ainv = svdinv(A)
    scale = sqrt(abs(diag(A)) + 1e-16);
    Anorm = diag(scale.^-1) * A * diag(scale.^-1);
                
    [U, S, V] = svd(Anorm);
    pos = diag(S) > max(diag(S)) / 1e6;

    Anorminv = V(:, pos) * diag(diag(S(pos, pos)).^-1) * U(:, pos)';
    Ainv = diag(scale.^-1) * Anorminv * diag(scale.^-1);
end