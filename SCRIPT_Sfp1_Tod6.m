% Setup
clearvars
close all

files = struct2cell(dir('Sfp1_Tod6'));                                      % all data sets
files = files(1, 3:end)';

t = 1:80;
Nm = 1000;

rng(0)


%% Fit flexible GPs and save results
for idx = 1:length(files)
    name = files{idx, 1}(1:end-11);
    disp(name)

    data = xlsread(sprintf('Sfp1_Tod6/%s', files{idx, 1}));                      % read data

    N = max(data(:, 11));                                                   % subsample, reshape, remove NaN, and normalize
    ratio = data(:, 7);
    data = reshape(ratio, 80, N);
    data = (data - mean(data)) ./ std(data);
    data = reshape(data(~isnan(data)), size(data, 1), []);
    data = data(t, 1:min(Nm, end));

    model = PoEmodel(t, data);                                              % define and fit GP model
    model.fit();

    files{idx, 2} = name;                                                   % organize results in cell array
    files{idx, 3} = data;
    files{idx, 4} = model;
end
                                                                            % convert to struct
separate = cell2struct(files, {'file', 'name', 'data', 'model'}, 2);
save('fitted2.mat', 'separate', 't');                                             % save for later use


%% Compare
close all
load('fitted2.mat');                                                         % load separate data and models

comparisons = [ 9*ones(8, 1)   (1:8)';                                      % specify desired comparisons
               19*ones(8, 1) (10:17)';                                      % currently: WT vs {mutants} for Tod6 and Sfp1 separately
                          10       16;                                      % + Tod6 6A vs Sch92d3e
                           5        4];                                     % + Sfp1 Pib2KO vs G1QL_G2SL_Pib2KO

combined2 = combined
combined = cell(length(comparisons), 4);

for idx = 1:4
    first = comparisons(idx, 1);                                            % first and second data counters
    second = comparisons(idx, 2);

    data1 = separate(first).data;
    data2 = separate(second).data;

    nseries = min(size(data1, 2), size(data2, 2));
    data1 = data1(:, 1:nseries);
    data2 = data2(:, 1:nseries);
    model = PoEmodel(t, [data1 data2]);
    model.fit();

    combined{idx, 1} = first;
    combined{idx, 2} = second;
    combined{idx, 3} = data1;
    combined{idx, 4} = data2;
    combined{idx, 5} = model;
    combined{idx, 6} = separate(first).model;
    combined{idx, 7} = separate(second).model;
end

combined = cell2struct(combined, {'first', 'second', 'data1', 'data2', 'model', 'model1', 'model2'}, 2);
save('fitted2.mat', 'separate', 'combined');                                 % save for later use


%% Posterior indicator expectations
close all
load('fitted2.mat');

[~, g9] = sgolay(3, 9);
t_filtered = [1:8:80 80];
diagonals = repmat(g9(:, 1)', 80, 1);
filter = spdiags(diagonals, -4:4, 80, 80);
filter = filter(t_filtered, :);
filter(1, 2:5) = filter(1, 2:5)*2;
filter(end, end-4:end-1) = filter(end, end-4:end-1)*2;
filter = full(filter);

% t_filtered = [1:8:80 80]
% filter = eye(80);
% filter = full(filter(t_filtered, :));

for idx = 1:length(combined)
    model1 = combined(idx).model1;
    model2 = combined(idx).model2;
    model = combined(idx).model;

    a = 1.01;
    b = 1.01;
    
    c = 2;
    d = 1.01;
    
    % gam_q_opt = marginal_likelihood(a, b, c, d, model, model1, model2);
    % gam_q_opt(3) = .26;
    gam_q_opt = [0 0 .25 .85];
    [E, MAP] = probabilities([.25 .85], model, model1, model2, filter);
    
    plot_filtered(model, model1, model2, E, filter);

    % tl = gcf().Children;
    % tl.Title.Interpreter = 'none';
    % tl.Title.String = strcat(separate([combined(idx).first]).name, " vs ", ...
    %                          separate([combined(idx).second]).name, ...
    %                          sprintf(" (γ_{opt} = [%.3f, %.3f], q_{opt} = %.3f", ...
    %                                  gam_q_opt(1), gam_q_opt(2), gam_q_opt(3)));

    file = strcat("figures/Sfp1_Tod6/", separate([combined(idx).first]).name, "_VS_", ...
                              separate([combined(idx).second]).name);
    exportgraphics(gcf, file, ContentType='vector')
end

compared = [num2cell(1:18)' {separate([combined.first]).name}' {separate([combined.second]).name}']