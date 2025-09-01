function data = generate(t, R, f, fA, fB, Z, hyperpar, smooth_function)
    model = PoEmodel(t, t);

    la = model.flexible_covariance(hyperpar.zeta);
    laA = model.flexible_covariance(hyperpar.zetaA);
    laB = model.flexible_covariance(hyperpar.zetaB);

    indep = diag(Z);
    shared = eye(length(Z)) - diag(Z);

    meanA = indep*fA' + shared*f';
    meanB = indep*fB' + shared*f';

    covA = indep*laA*indep + shared*la*shared;
    covB = indep*laB*indep + shared*la*shared;

    data = struct();
    data.YA = mvnrnd(meanA', covA, R);
    data.YB = mvnrnd(meanB', covB, R);
end