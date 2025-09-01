% Synthetic 1
% clearvars
% close all

zeta = [-3*ones(1, 7) 1.45*ones(1, 7) -10];
hyperpar = struct('zeta', zeta, 'zetaA', zeta, 'zetaB', zeta);

fA = zeros(1, 11);
f = fA;
fB = [zeros(1, 5), linspace(0, .125, 6)];

t = 0:8:80;

Z = [1 1 1 1 1 1 1 1 1 1 1];
Ztrue = [zeros(1, 6) ones(1, 5)];

bspline = BSpline(4, [.25 .5 .75], (t - t(1))/range(t));


R_values = [5 15 50];
seeds = 1:10;
expectations = zeros(length(t), length(seeds), length(R_values));
aurocs = zeros(length(seeds), length(R_values));

% for Ridx = 1:length(R_values)%, for seed = seeds
%     seeds = [4 7 9];
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
%     save('simulation/meangradual.mat')
% 
%     plot_new(model, model1, model2, E, ["fit" "mprior" "mpost"]);
% 
%     % tl = gcf().Children;
%     % tl.Title.Interpreter = 'TeX';
%     % % title(tl, sprintf("γ_{opt} = [%.3f, %.3f], q_{opt} = %.3f", ...
%     % %                   gam_q_opt(1), gam_q_opt(2), gam_q_opt(3)))
%     % title(tl, "{\itIn silico} data: \bfMean difference")
% 
%     file = sprintf("figures/insilico_meangradual_R%d_%d.pdf", R, seed);
% 
%     exportgraphics(gcf, file, ContentType='vector')
% end%, end


%%
% Ztrue = [zeros(1, 12) ones(1, 4)];
% for Ridx = 1:3, for seed = seeds 
%     E = expectations(:, seed, Ridx)';
%     aurocs2(seed, Ridx) = AUCROC(E, [zeros(1, 12) ones(1, 4)]);
% end, end
% 
% aurocs
% aurocs2

for Ridx = 1:length(R_values)%, for seed = seeds
    seeds = [4 7 9];
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
    % save('simulation/meangradual.mat')

    plot_new(model, model1, model2, E, [fA; fB/scale]);
    
    file = sprintf("figures/insilico_meangradual_R%d_%d.pdf", R, seed);
    
    exportgraphics(gcf, file, ContentType='vector')
end%, end
