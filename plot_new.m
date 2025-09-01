function plot_new(model, modelA, modelB, E, fAB, cAB)
    if nargin < 5, fAB = []; cAB = []; end
    if nargin < 6, cAB = []; end
    % figure(Position=[100 100 1600 1000])
    % tiledlayout(2, 3+(exist('E', 'var') && ~isempty(E)))

    figure(Position=[100 100 950 350])
    tl = tiledlayout(1, 4, TileSpacing='compact');
    tl.XLabel.String = "{\itt} (relative time index)";
    tl.YLabel.String = "\itY";
    
    plot_data(modelA, modelB, fAB, cAB)
    plot_fit(model, modelA, modelB)
    plot_mpost(model, modelA, modelB, fAB, cAB)

    plot_indicators(model, modelA, modelB, E)
end


function plot_data(modelA, modelB, fAB, cAB)
    t = modelA.t;
    data1 = modelA.y;
    data2 = modelB.y;
    R1 = modelA.R;
    R2 = modelB.R;

    nexttile(1)
    hold off
    plot(t, data1, 'LineWidth', 3, 'Color', [0 0.4470 0.7410 .5/sqrt(R1)])
    hold on
    if ~isempty(fAB)
        plot(t, fAB(1, :), 'k', LineWidth=3, Color=[0 0 0 .5]);
    end
    if ~isempty(cAB)
        plot(linspace(0, 80, 81), -cAB(1, :), 'k:', LineWidth=3, Color=[0 0 0 .5]);
        plot(linspace(0, 80, 81), cAB(1, :), 'k:', LineWidth=3, Color=[0 0 0 .5]);
    end
    
    ylim([-4 4])
    lim1 = ylim;
    % xlabel('\itt')
    % ylabel('\itY')
    ax = gca;
    ax.Color = [.975 .975 .975];
    % set(ax, 'XTickLabel', [])
    xlim([0 t(end)])
    xticks(t(end) .* (0:.25:1))
    xticklabels(["0" "0.25" "0.5" "0.75" "1"])
    xtickangle(0)
    
    nexttile(2)
    hold off
    plot(t, data2, 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980 .5/sqrt(R2)])
    hold on
    if ~isempty(fAB)
        plot(t, fAB(2, :), 'k', LineWidth=3, Color=[0 0 0 .5]);
    end
    if ~isempty(cAB)
        plot(linspace(0, 80, 81), -cAB(2, :), 'k:', LineWidth=3, Color=[0 0 0 .5]);
        plot(linspace(0, 80, 81), cAB(2, :), 'k:', LineWidth=3, Color=[0 0 0 .5]);
    end
    
    ylim([-4 4])
    lim2 = ylim;
    % xlabel('\itt')
    % ylabel('\itY')
    ax = gca;
    ax.Color = [.975 .975 .975];
    % set(ax, 'XTickLabel', [])
    set(ax, 'YTickLabel', [])
    xlim([0 t(end)])
    xticks(t(end) .* (0:.25:1))
    xticklabels(["0" "0.25" "0.5" "0.75" "1"])
    
    nexttile(3)
    hold off
    plot(t, data1, 'LineWidth', 4, 'Color', [0 0.4470 0.7410 .5/sqrt(R1)/2])
    hold on
                    ylim([-4 4])
    plot(t, data2, 'LineWidth', 4, 'Color', [0.8500 0.3250 0.0980 .5/sqrt(R2)/2])
    lim3 = ylim;
    % xlabel('\itt')
    % ylabel('\itY')
    ax = gca;
    ax.Color = [.975 .975 .975];
    xlim([0 t(end)])
    set(ax, 'YTickLabel', [])
    % xticks([1 xticks])
    xticks(t(end) .* (0:.25:1))
    xticklabels(["0" "0.25" "0.5" "0.75" "1"])
    
    lim = max(max(abs([lim1; lim2; lim3]))) * [-1 1] + 0;
    for idx = 1:3
        nexttile(idx)
        ylim(lim)
    end
end


function plot_fit(model, modelA, modelB)
    t = model.t;
    t_fine = linspace(t(1), t(end), 159);
    c = 2;
    
    mu1  = interp1(t, modelA.muf, t_fine)';
    mu2  = interp1(t, modelB.muf, t_fine)';
    mu  = interp1(t, model.muf, t_fine)';
    
    s1 = sqrt(diag(modelA.flexible_covariance(modelA.ze, t_fine)));
    s2 = sqrt(diag(modelB.flexible_covariance(modelB.ze, t_fine)));
    s = sqrt(diag(model.flexible_covariance(model.ze, t_fine)));
    
    nexttile(1)                                               % add means +- sds to plots
    plot(t_fine, [mu1, mu1 + c*s1, mu1 - c*s1], ':', 'LineWidth', 3, 'Color', [0 0.4470 0.7410 .5])
        
    nexttile(2)
    plot(t_fine, [mu2, mu2 + c*s2, mu2 - c*s2], ':', 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980 .5])
    
    nexttile(3)
    plot(t_fine, [mu, mu + c*s, mu - c*s], ':', 'LineWidth', 3, 'Color', [0.9290 0.6940 0.1250 .75])
end


function plot_mpost(model, modelA, modelB, fAB, cAB)
    t = model.t;
    c = 2;

    mu1 = modelA.muf;
    s1  = modelA.sf;
    
    mu2 = modelB.muf;
    s2  = modelB.sf;
    
    mu = model.muf;
    s  = model.sf;
    
    nexttile(1)                                               % add means +- sds to plots
    plot(t, [mu1, mu1 + c*s1, mu1 - c*s1], '-', 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410 .75])
    plot(t, 0*t+100, 'LineWidth', 3, 'Color', [0 0.4470 0.7410 .25])
    ax = gca;
    h1 = ax.Children(1);
    h2 = ax.Children(5);
    h3 = ax.Children(2);
    h4 = ax.Children(8);

    if ~isempty(fAB)
        legend([h4, h1, h2, h3], 'Ground truth mean', 'Data', 'Measurement GP', 'Latent GP posterior', location='north')
    end
    if ~isempty(cAB)
        legend([h4, h1, h2, h3], 'Ground truth cov.', 'Data', 'Measurement GP', 'Latent GP posterior', location='north')
    end

    nexttile(2)
    plot(t, [mu2, mu2 + c*s2, mu2 - c*s2], '-', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980 .75])
    
    nexttile(3)
    plot(t, [mu, mu + c*s, mu - c*s], '-', 'LineWidth', 1.5, 'Color', [0.9290 0.6940 0.1250 .875])
end


function plot_indicators(model, modelA, modelB, E)
    % col_EZ = 'w';
    col_EZ = 'k';

    t = model.t;
    t_fine = linspace(t(1), t(end), 159);
    data1 = modelA.y;
    data2 = modelB.y;
    c = 2;

    timepoint = t';
    p_Z = E';
    table(timepoint, p_Z)
    
    


    mu1  = interp1(t, modelA.muf, t_fine)';
    mu2  = interp1(t, modelB.muf, t_fine)';
    mu  = interp1(t, model.muf, t_fine)';
    
    s1 = sqrt(diag(modelA.flexible_covariance(modelA.ze, t_fine)));
    s2 = sqrt(diag(modelB.flexible_covariance(modelB.ze, t_fine)));
    s = sqrt(diag(model.flexible_covariance(model.ze, t_fine)));

    nexttile(1)
    lim = ylim;
    tile = nexttile(4);

    colors = summer(256);
    colors = 2/3 + colors/3;
    colors = colors(end:-1:1, :);
    
    yyaxis left
    hold on
    for n = 1:length(t)
        col = colors(round(1 + p_Z(n) * 255), :);
        xl = t(n) - (t(n) - t(max(1, n-1)))/2;
        xr = t(n) + (t(min(end, n+1)) - t(n))/2;
        patch([xl xr xr xl], [lim(1) lim(1) lim(2) lim(2)], col, FaceAlpha=1, EdgeAlpha=0)
    end

    % patch([1 1 1 1], [lim(1) lim(1) lim(2) lim(2)], colors(1, :), FaceAlpha=1, EdgeColor='k')
    % patch([1 t(end) t(end) 1], [lim(1) lim(1) lim(1) lim(1)], colors(end, :), FaceAlpha=1, EdgeColor='k')
    set(gca, 'Layer', 'top')
    
    % plot(t, data1, '-', 'LineWidth', 4, 'Color', [0 0.4470 0.7410 .05], 'MarkerSize', eps)
    % plot(t, data2, '-', 'LineWidth', 4, 'Color', [0.8500 0.3250 0.0980 .05], 'MarkerSize', eps)
    plot(t_fine, [mu1, mu1 + c*s1, mu1 - c*s1], ':', 'LineWidth', 3, 'Color', [0 0.4470 0.7410 .5])
    plot(t_fine, [mu2, mu2 + c*s2, mu2 - c*s2], ':', 'LineWidth', 3, 'Color', [0.8500 0.3250 0.0980 .5])
    plot(t_fine, [mu, mu + c*s, mu - c*s], ':', 'LineWidth', 3, 'Color', [0.9290 0.6940 0.1250 .5])
    % plot(t, [mu + s, mu - s], ':', 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250 .6])
    
    xlim([0 t(end)])
    xticks(t(end) .* (0:.25:1))
    xticklabels(["0" "0.25" "0.5" "0.75" "1"])
    % xticks([1 xticks])
    ylim(lim)
    % xlabel('\itt')
    % ylabel('\itY')
    set(gca, 'YTickLabel', [])

    yyaxis right
    tile.YAxis(1).Color = 'k';
    tile.YAxis(2).Color = "#A2142F";
    logit_Z = log(p_Z ./ max(eps, (1 - p_Z)));
    plot(t, logit_Z, Color="#A2142F", LineWidth=2)
    lim = max(abs(ylim));
    ylim(lim*[-1 1])
    ylabel('logit p(\itZ_n=1)')
    % ax = gca;
    % ax.YLabel.Position(1) = t(end)*1.15;

    % yline((lim(2)-1)/2 * [-1 1] + 1)
    % plot(timepoint, (lim(2) - 1) * (p_Z - 1/2) + 1, [col_EZ '--'], 'LineWidth', 2)
    % plot([20 30], 1.1*(lim(2)-1)/2 * [1 1] + 1, [col_EZ '--'], 'LineWidth', 2)
    % text(19, 1.1*(lim(2)-1)/2 + 1, 'p(Zn)=1', 'HorizontalAlignment', 'right')
    % text(19, -1.1*(lim(2)-1)/2 + 1, 'p(Zn)=0', 'HorizontalAlignment', 'right')
    % 
    % logit_Z = log(p_Z ./ max(eps, (1 - p_Z)));
    % max_logit = ceil(max(abs(logit_Z)));
    % plot(timepoint, (lim(2)-1)/2 * (logit_Z/max_logit) + 1, [col_EZ ':'], 'LineWidth', 2)
    % text(60, 1.1*(lim(2)-1)/2 + 1, sprintf('logit(Zn)=%d', max_logit))
    % text(60, -1.1*(lim(2)-1)/2 + 1, sprintf('logit(Zn)=%d', -max_logit))
    % plot([50 59], 1.1*(lim(2)-1)/2 * [1 1] + 1, [col_EZ ':'], 'LineWidth', 2)
    
    colormap(tile, colors)
    cb = colorbar(tile);
    cb.Ticks = [0 1];
    cb.TickLabels = {'0', '1'};
    cb.Label.String = 'p({\itZ}_n=1)';
    cb.Label.Position(1) = 1.2;
    disp(1)
end