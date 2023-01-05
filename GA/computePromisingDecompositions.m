clear;
close all;
% clc;

%% 

restoredefaultpath();
system_name = 'biped2d';
addpath('utils');
addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
addpath('../iLQG_boxQP/iLQG utilities/decomposition_count');
% load(strcat('../iLQG_boxQP/new_systems/', system_name, '/', system_name, 'System.mat'));
load(strcat('../MCTS/data/', system_name, 'System.mat'));
% sys.I(2,2) = 0.5*(sys.I(1,1) + sys.I(3,3));
% sys.gamma_ = 1;
sys.num_action_samples = [5, 5, 2, 2];
sys.max_iter = 2000;
sys.state_bounds(1,:) = [0.92, 1];
sys.state_bounds(2,:) = [0.2, 0.3] + pi/2;

% Add function to compute linearized dynamics
sys.fxfu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
fxfu = sys.fxfu_func(sys.l_point, sys.u0);
sys.A = fxfu(:, 1:sys.X_DIMS);
sys.B = fxfu(:, (1+sys.X_DIMS):end);

% Add function to compute LQR metric
sys.lambda_ = (1 - sys.gamma_) / sys.dt;
[K_j, S_j, ~] = lqr(sys.A - (sys.lambda_ / 2) * eye(size(sys.A,1)), sys.B, sys.Q, sys.R, zeros(size(sys.B)));
% sys.S = sym('S', [sys.X_DIMS, sys.X_DIMS]);
% sys.S = tril(sys.S,0) + tril(sys.S,-1).';
% sys.a = sym('a', [sys.X_DIMS, 1]);
% sys.b = sym('b', [sys.X_DIMS, 1]);
% err_lqr = sys.x.'*sys.S*sys.x - sys.x.'*S_j*sys.x;
% for ii=1:1:sys.X_DIMS
%     err_lqr = int(err_lqr, sys.x(ii), [sys.a(ii), sys.b(ii)]);
% end
% sys.err_lqr_func = matlabFunction(err_lqr, 'Vars', {sys.S, sys.a, sys.b});
sys.err_lqr_func = @(S,a,b) definite_integral_parabola(S-S_j,a,b);
% sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1), 1);
sys.da = definite_integral_parabola(S_j, sys.state_bounds(:,1)-sys.l_point, sys.state_bounds(:,2)-sys.l_point);

sys.measure_func = @(err_lqr, err_compute) (1 - exp(-err_lqr)) .* err_compute;
% sys.measure_func = @(err_lqr, err_compute) (min(20.0, err_lqr) / 20.0) .* err_compute;

% sys.X_DIMS_MASKED = [eye(sys.X_DIMS/2), eye(sys.X_DIMS/2)];
sys.X_DIMS_MASKED = eye(sys.X_DIMS);
% sys.num_points = sys.numPoints;
% sys.num_action_samples = 10*ones(1, sys.U_DIMS);
% sys.max_iter = 2000;
% sys.max_policy_iter = 100;
% sys.num_action_samples = 40*ones(1, sys.U_DIMS);
% sys.max_iter = 5000;
% sys.max_policy_iter = 800;
% sys.num_points = [7,7,7,35,7,7,7,11,11,35];

nvars = 2*sys.U_DIMS + sys.X_DIMS*sys.U_DIMS;
lb = [zeros(1, sys.U_DIMS), ones(1, sys.U_DIMS), zeros(1, sys.U_DIMS*sys.X_DIMS)];
ub = [sys.U_DIMS*ones(1, sys.U_DIMS*2), ones(1, sys.U_DIMS*sys.X_DIMS)];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = @(x) constraints(reshape(x(1:2*sys.U_DIMS), sys.U_DIMS, 2), ...
                           reshape(x((2*sys.U_DIMS + 1):(2*sys.U_DIMS + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS));
IntCon = [];
options.Display = 'iter';
options.PopulationSize = 200;
options.CrossoverFraction = 0.0;
options.EliteCount = 0.2*options.PopulationSize;

options.CreationFcn = @(nvars, fitness_fcn, options) ...
                        uniform_input_tree_sampling(sys, options.PopulationSize);
options.CrossoverFcn = @(parents, options, nvars, fitness_fcn, unused, population) ...
                         crossoverfunction(sys, parents, options, nvars, fitness_fcn, unused, population);
% options.MutationFcn = @(parents, options, nvars, fitness_fcn, state, score, population) ...
%                          mutationfunction(sys, parents, options, nvars, fitness_fcn, state, score, population);
options.MutationFcn = @(parents, options, nvars, fitness_fcn, state, score, population) ...
                         mutationfunctionstatemasked(sys, parents, options, nvars, fitness_fcn, state, score, population);
MaxGATime = 150;
MaxTotalTime = 150;
num_runs = 5;
ga_stats = zeros(num_runs, 7);
ga_decompositions = cell(num_runs, 2);
decoded_ga_decompositions = cell(num_runs, 3); % action_coupling, action_dependence, state_dependence

for rr=1:1:num_runs
    
    sys.decompositionlist = containers.Map();
    fun = @(x) computeJointMeasure(sys, reshape(x(1:2*sys.U_DIMS), sys.U_DIMS, 2), ...
                                   reshape(x((2*sys.U_DIMS + 1):(2*sys.U_DIMS + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS));
    ga_solutions = cell(0, 6); % x, err_joint, exitflag, output, population, scores
    tic;
    while((MaxTotalTime - toc) > 0)
        options.InitialPopulation = uniform_input_tree_sampling(sys, options.PopulationSize);
        options.MaxTime = min(MaxGATime, MaxTotalTime - toc);
        [ga_solutions{(end+1), :}] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    end

    % Extract decomposition
    num_to_extract = 10;
    final_population = cell2mat(ga_solutions(:, end-1));
    final_scores     = cell2mat(ga_solutions(:, end));

    [~, unique_population_id] = extract_decompositions_from_population(sys, final_population);
    final_population = final_population(unique_population_id, :);
    final_scores     = final_scores(unique_population_id, :);
    [~, final_scores_order] = sort(final_scores);

    final_population = final_population(final_scores_order(1:min(num_to_extract, size(final_population, 1))), :);
    final_scores     = final_scores(final_scores_order(1:min(num_to_extract, size(final_scores, 1))), :);

    [cumm_population_encoding, cumm_population_ids] = unique(cell2mat(sys.decompositionlist.keys'), ...
                                                             'rows', 'stable');
    cumm_population_measure = cell2mat(sys.decompositionlist.values');
    cumm_population_measure = cumm_population_measure(cumm_population_ids,:);

    [~, best_population_ids] = sort(cumm_population_measure(:,1));
    best_population_ids = best_population_ids(1:num_to_extract);
    best_population_measure = cumm_population_measure(best_population_ids, :);
    best_population = cumm_population_encoding(best_population_ids, :);

    best_measure = min(best_population_measure(:,1));
    mean_measure = mean(best_population_measure(:,1));
    std_measure = std(best_population_measure(:,1));

    best_lqr_measure = min(best_population_measure(:,2));
    mean_lqr_measure = mean(best_population_measure(:,2));
    std_lqr_measure = std(best_population_measure(:,2));
    
    ga_stats(rr, 1) = best_measure;
    ga_stats(rr, 2) = mean_measure;
    ga_stats(rr, 3) = std_measure;
    
    ga_stats(rr, 4) = best_lqr_measure;
    ga_stats(rr, 5) = mean_lqr_measure;
    ga_stats(rr, 6) = std_lqr_measure;
    
    ga_stats(rr, 7) = size(sys.decompositionlist, 1);
    
    ga_decompositions{rr, 1} = best_population;
    ga_decompositions{rr, 2} = best_population_measure;
    
    disp('GA');
    disp(strcat('Best Measure : ', num2str(best_measure), ...
                ', Mean Measure : ', num2str(mean_measure), ...
                ', Std Measure : ', num2str(std_measure)));
    disp(strcat('Best LQR Measure : ', num2str(best_lqr_measure), ...
                ', Mean LQR Measure : ', num2str(mean_lqr_measure), ...
                ', Std LQR Measure : ', num2str(std_lqr_measure)));
    disp(strcat('', 'Number of nodes : ', num2str(size(sys.decompositionlist, 1))));
    
    % Decode decompositions
    decoded_ga_decompositions{rr, 1} = zeros(sys.U_DIMS, sys.U_DIMS, num_to_extract);
    decoded_ga_decompositions{rr, 2} = zeros(sys.U_DIMS, sys.U_DIMS, num_to_extract);
    decoded_ga_decompositions{rr, 3} = zeros(sys.U_DIMS, sys.X_DIMS, num_to_extract);
    for dd=1:1:num_to_extract
        decoded_ga_decompositions{rr, 1}(:,:,dd) ...
         = reshape(ga_decompositions{rr,1}(dd, 1:sys.U_DIMS^2) - 48, sys.U_DIMS, sys.U_DIMS);
        decoded_ga_decompositions{rr, 2}(:,:,dd) ...
         = reshape(ga_decompositions{rr,1}(dd, (1+sys.U_DIMS^2):(2*sys.U_DIMS^2)) - 48, sys.U_DIMS, sys.U_DIMS);
        decoded_ga_decompositions{rr, 3}(:,:,dd) ...
         = reshape(ga_decompositions{rr,1}(dd, (1+2*sys.U_DIMS^2):end) - 48, sys.U_DIMS, sys.X_DIMS);
    end
end

%% Random Sampling

num_to_sample = 1000;
random_stats = zeros(num_runs, 7);
random_decompositions = cell(num_runs, 3);
for rr=1:1:num_runs
    sys.decompositionlist = containers.Map();
    sampled_population = cell(0, 2);
    tic;
    while ((MaxTotalTime - toc) > 0)

        new_population = uniform_input_tree_sampling(sys, num_to_sample);
        new_population = num2cell(new_population, 2);
        new_population_measure = cellfun(@(x) computeJointMeasure(sys, ...
                                                                  reshape(x(1:(2*sys.U_DIMS)), sys.U_DIMS, 2), ...
                                                                  reshape(x((1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS)), ...
                                         new_population, 'UniformOutput', false);
        sampled_population = cat(1, sampled_population, cat(2, new_population, new_population_measure));
        if (size(sampled_population, 1) > 1e4)
            [~, unique_sampled_population_id] = extract_decompositions_from_population(sys, cell2mat(sampled_population(:,1)));
            sampled_population = sampled_population(unique_sampled_population_id, :);
        end
    end
    
    num_to_extract = 10;
    [~, unique_sampled_population_id] = extract_decompositions_from_population(sys, cell2mat(sampled_population(:,1)));
    sampled_population_measure = cell2mat(sampled_population(unique_sampled_population_id, 2));
    sampled_population = cell2mat(sampled_population(unique_sampled_population_id, 1));

    [sampled_population_measure, sampled_population_order] = sort(sampled_population_measure);
    sampled_population = sampled_population(sampled_population_order(1:min(num_to_extract, size(sampled_population, 1))), :);
    sampled_population_measure = sampled_population_measure(1:size(sampled_population, 1), :);

    sampled_population_err_lqr = zeros(num_to_extract, 1);
    for dd=1:1:num_to_extract
        sampled_population_err_lqr(dd) = computeLQRMeasure(sys, ...
                                                           reshape(sampled_population(dd, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2), ...
                                                           reshape(sampled_population(dd, (1+(2*sys.U_DIMS)):end), sys.U_DIMS, sys.X_DIMS));
    end
    
    random_stats(rr, 1) = min(sampled_population_measure);
    random_stats(rr, 2) = mean(sampled_population_measure);
    random_stats(rr, 3) = std(sampled_population_measure);
    
    random_stats(rr, 4) = min(sampled_population_err_lqr);
    random_stats(rr, 5) = mean(sampled_population_err_lqr);
    random_stats(rr, 6) = std(sampled_population_err_lqr);
    
    random_stats(rr, 7) = size(sys.decompositionlist, 1);
    
    random_decompositions{rr, 1} = sampled_population;
    random_decompositions{rr, 2} = sampled_population_measure;
    random_decompositions{rr, 3} = sampled_population_err_lqr;
    
    disp('Random');
    disp(strcat('Best Measure : ', num2str(min(sampled_population_measure)), ...
                ', Mean Measure : ', num2str(mean(sampled_population_measure)), ...
                ', Std Measure : ', num2str(std(sampled_population_measure))));
    disp(strcat('Best LQR Measure : ', num2str(min(sampled_population_err_lqr)), ...
                ', Mean LQR Measure : ', num2str(mean(sampled_population_err_lqr)), ...
                ', Std LQR Measure : ', num2str(std(sampled_population_err_lqr))));
    disp(strcat('', 'Number of nodes : ', num2str(size(sys.decompositionlist, 1))));
end
sys.decompositionlist = [];

% save(strcat('data/', system_name,'_GA_NoCrossover_RandomSampled_',num2str(MaxTotalTime),'_run5.mat'), ...
%      'sys', 'final_population', 'final_scores', 'best_population', 'best_population_measure', 'sampled_population', 'sampled_population_measure', 'sampled_population_err_lqr', 'MaxTotalTime', 'MaxGATime', 'num_to_extract');

save(strcat('data/', system_name,'_GA_NoCrossover_UniformSampled_explqrobj',num2str(MaxTotalTime),'_newactionsample.mat'), ...
     'sys', 'ga_stats', 'ga_decompositions', 'decoded_ga_decompositions', 'random_stats', 'random_decompositions', 'MaxTotalTime', 'MaxGATime', 'num_to_extract');
