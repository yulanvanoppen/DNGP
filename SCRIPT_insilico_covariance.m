% Synthetic 1
clearvars
close all

t = 0:8:80;
model0 = PoEmodel(t, 0*t);

zeta = [-3*ones(1, 7) 1.5*ones(1, 7) -10];
zeta2 = zeta;

X = model0.bspline.cubic_basis_base;

amplitude1_aug = exp(X * zeta(1:7)') .* (1 + .75 * [zeros(1, 3) ones(1, 3), zeros(1, 5)]');
lamp1 = log(amplitude1_aug);
zeta2(1:7) = (X' * X) \ (X' * lamp1);

lscale1_aug = exp(X * zeta(8:14)') .* (1 - .75 * [zeros(1, 3) ones(1, 3), zeros(1, 5)]');
llscale1 = log(lscale1_aug);
zeta2(8:14) = (X' * X) \ (X' * llscale1);

hyperpar = struct('zeta', zeta, 'zetaA', zeta, 'zetaB', zeta2);

fA = zeros(1, 11);
f = fA;
fB = fA;

Z = [1 1 1 1 1 1 1 1 1 1 1];
Ztrue = [zeros(1, 3) ones(1, 3), zeros(1, 5)];

bspline = BSpline(4, [.25 .5 .75], (t - t(1))/range(t));


R_values = [5 15 50];
seeds = 1:10;
expectations = zeros(length(t), length(seeds), length(R_values));
aurocs = zeros(length(seeds), length(R_values));

% for Ridx = 1:length(R_values)%, for seed = seeds
%     seeds = [5 8 5];
%     seed = seeds(Ridx);
% 
%     rng(seed)
%     R = R_values(Ridx);
% 
% 
%     %% Generate data
%     data = generate(t, R, f, fA, fB, Z, hyperpar, bspline);
%     data1 = data.YA';
%     data2 = data.YB';
% 
%     scale = std([data1 data2], 0, 'all');
%     data1 = data1 ./ scale;
%     data2 = data2 ./ scale;
% 
% 
%     %% Fit hyperparameters
%     model1 = PoEmodel(t, data1);
%     model1.fit(50, 1e-3);
%     out1 = struct('theta', model1.theta_opt(end, :), 'zeta', model1.zeta_opt(end, :), ...
%                      'muf', model1.muf, 'Sf', model1.Sf, 'sf', model1.sf);
% 
%     model2 = PoEmodel(t, data2);
%     model2.fit(50, 1e-3);
%     out2 = struct('theta', model2.theta_opt(end, :), 'zeta', model2.zeta_opt(end, :), ...
%                      'muf', model2.muf, 'Sf', model2.Sf, 'sf', model2.sf);
% 
%     model = PoEmodel(t, [data1 data2]);
%     model.fit(50, 1e-3);
%     out = struct('theta', model.theta_opt(end, :), 'zeta', model.zeta_opt(end, :), ...
%                      'muf', model.muf, 'Sf', model.Sf, 'sf', model.sf);
% 
% 
%     %% Posterior probabilities and plot
%     [E, MAP] = probabilities([.25 .85], model, model1, model2);
%     aurocs(seed, Ridx) = AUCROC(E, Ztrue);
%     expectations(:, seed, Ridx) = E;
% 
%     close all
%     save('simulation/covariance.mat')
% 
%     plot_new(model, model1, model2, E, ["fit" "mprior" "mpost"]);
% 
%     % tl = gcf().Children;
%     % tl.Title.Interpreter = 'TeX';
%     % % title(tl, sprintf("γ_{opt} = [%.3f, %.3f], q_{opt} = %.3f", ...
%     % %                   gam_q_opt(1), gam_q_opt(2), gam_q_opt(3)))
%     % title(tl, "{\itIn silico} data: \bfCovariance difference")
% 
%     file = sprintf("figures/insilico_covariance_R%d_%d.pdf", R, seed);
% 
%     exportgraphics(gcf, file, ContentType='vector')
% end%, end



%%
cov1 = interp1(0:8:80, exp(X * zeta(1:7)')', 0:80);
cov2 = interp1(0:8:80, amplitude1_aug', 0:80);
cov2 = conv(cov2, ones(1, 7)/7, 'same');
cov2([1:7 end-6:end]) = cov1([1:7 end-6:end]);

for Ridx = 1:length(R_values)%, for seed = seeds
    seeds = [5 8 5];
    seed = seeds(Ridx);

    rng(seed)
    R = R_values(Ridx);


    %% Generate data
    data = generate(t, R, f, fA, fB, Z, hyperpar, bspline);
    data1 = data.YA';
    data2 = data.YB';
    
    scale = std([data1 data2], 0, 'all');
    data1 = data1 ./ scale;
    data2 = data2 ./ scale;
    
    
    %% Fit hyperparameters
    model1 = PoEmodel(t, data1);
    model1.fit(50, 1e-3);
    out1 = struct('theta', model1.theta_opt(end, :), 'zeta', model1.zeta_opt(end, :), ...
                     'muf', model1.muf, 'Sf', model1.Sf, 'sf', model1.sf);
    
    model2 = PoEmodel(t, data2);
    model2.fit(50, 1e-3);
    out2 = struct('theta', model2.theta_opt(end, :), 'zeta', model2.zeta_opt(end, :), ...
                     'muf', model2.muf, 'Sf', model2.Sf, 'sf', model2.sf);
    
    model = PoEmodel(t, [data1 data2]);
    model.fit(50, 1e-3);
    out = struct('theta', model.theta_opt(end, :), 'zeta', model.zeta_opt(end, :), ...
                     'muf', model.muf, 'Sf', model.Sf, 'sf', model.sf);
    
    
    %% Posterior probabilities and plot
    [E, MAP] = probabilities([.25 .85], model, model1, model2);
    
    close all
    % save('simulation/covariance.mat')

    plot_new(model, model1, model2, E, [], 2/scale*[cov1; cov2]);
    
    file = sprintf("figures/insilico_covariance_R%d_%d.pdf", R, seed);
    
    exportgraphics(gcf, file, ContentType='vector')
end%, end
