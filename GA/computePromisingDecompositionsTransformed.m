clear;
close all;
% clc;

%% 

restoredefaultpath();
system_name = 'transformed_quadcopter_wload_reduced';
if (strcmp(system_name, 'quadcopter_double_transformed'))
    load_system_name = 'quadcopter_transformed';
elseif (strcmp(system_name, 'transformed_quadcopter'))
    load_system_name = 'quadcopter';
elseif (strcmp(system_name, 'transformed_quadcopter_wload'))
    load_system_name = 'quadcopter_wload';
elseif (strcmp(system_name, 'transformed_quadcopter_wload_reduced'))
    load_system_name = 'quadcopter_wload_reduced';
else
    load_system_name = split(system_name, '_');
    assert(length(load_system_name)>=2 && strcmp(load_system_name{end}, 'transformed'));
    load_system_name = split(system_name, '_transformed');
    load_system_name = load_system_name{1};
end
addpath('utils');
if (strcmp(system_name, 'unicycle_wflywheel_tilted_transformed'))
    addpath(strcat('../iLQG_boxQP/new_systems/', 'unicycle_wflywheel_transformed'));
else
    addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
end
addpath('../iLQG_boxQP/iLQG utilities/decomposition_count');
addpath('../../OptM');
% load(strcat('../MCTS/data/', load_system_name, 'System.mat'));
load(strcat('../iLQG_boxQP/new_systems/', load_system_name, '/', load_system_name, 'System.mat'));

% Compute transformation
% sys.I(2,2) = 0.5*(sys.I(1,1) + sys.I(3,3));
sys.Tx = eye(sys.X_DIMS);
sys.Tu = eye(sys.U_DIMS);
A = dynx(sys, zeros(sys.X_DIMS, 1), zeros(sys.U_DIMS, 1));
B = dynu(sys, zeros(sys.X_DIMS, 1), zeros(sys.U_DIMS, 1));
Q = sys.Q;
R = sys.R;
sys.lambda_ = (1 - sys.gamma_) / sys.dt;
[UK,SK,VK] = optimize_transform(A,B,Q,R,sys.lambda_);
sys.Tx = VK;
sys.Ty = sys.Tx';
sys.Tu = UK;
sys.Tv = sys.Tu';
sys.Q = sys.Tx' * sys.Q * sys.Tx;
sys.R = sys.Tu' * sys.R * sys.Tu;

% Add function to compute linearized dynamics
sys.fxfu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
fxfu = sys.fxfu_func(zeros(sys.X_DIMS, 1), zeros(sys.U_DIMS, 1));
sys.A = (sys.Tx') * A * sys.Tx;
sys.B = (sys.Tx') * B * sys.Tu;

% Add function to compute LQR metric
sys.lambda_ = (1 - sys.gamma_) / sys.dt;
[K_j, S_j, ~] = lqr(sys.A - (sys.lambda_ / 2) * eye(size(sys.A,1)), sys.B, sys.Q, sys.R, zeros(size(sys.B)));
[~, K_i, ~, info] = icare(sys.A - sys.lambda_/2*eye(sys.X_DIMS), sys.B, sys.Q, sys.R, [], [], []);
sys.err_lqr_func = @(S,a,b) definite_integral_parabola(sys.Tx*(S-S_j)*sys.Tx',a,b);
% sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1), 1);
sys.da = definite_integral_parabola(sys.Tx*S_j*sys.Tx', sys.state_bounds(:,1)-sys.l_point, sys.state_bounds(:,2)-sys.l_point);

sys.measure_func = @(err_lqr, err_compute) (1 - exp(-err_lqr)) .* err_compute;

sys.X_DIMS_MASKED = eye(sys.X_DIMS);

%%
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
MaxGATime = 1500;
MaxTotalTime = 1500;
num_runs = 1;
ga_stats = zeros(num_runs, 7);
ga_decompositions = cell(num_runs, 2);
decoded_ga_decompositions = cell(num_runs, 3); % action_coupling, action_dependence, state_dependence

for rr=1:1:num_runs
    
    sys.decompositionlist = containers.Map();
    fun = @(x) computeJointMeasureTransformed(sys, reshape(x(1:2*sys.U_DIMS), sys.U_DIMS, 2), ...
                                   reshape(x((2*sys.U_DIMS + 1):(2*sys.U_DIMS + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS));
    ga_solutions = cell(0, 6); % x, err_joint, exitflag, output, population, scores
    tic;
    while((MaxTotalTime - toc) > 0)
        options.InitialPopulation = uniform_input_tree_sampling(sys, options.PopulationSize);
        options.MaxTime = min(MaxGATime, MaxTotalTime - toc);
        [ga_solutions{(end+1), :}] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
    end

    % Extract decomposition
    num_to_extract = 100;
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
        new_population_measure = cellfun(@(x) computeJointMeasureTransformed(sys, ...
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

save(strcat('data/', system_name,'_GA_NoCrossover_UniformSampled_explqrobj',num2str(MaxTotalTime),'_optimizedTransform_rerun.mat'), ...
     'sys', 'ga_stats', 'ga_decompositions', 'decoded_ga_decompositions', 'random_stats', 'random_decompositions', 'MaxTotalTime', 'MaxGATime', 'num_to_extract');

function [UK,SK,VKopt] = optimize_transform(A,B,Q,R,lambda_)

[~,K,~,info] = icare(A - lambda_/2*eye(size(A,1)), B, Q, R, [], [], []);
[UK,SK,VK] = svd(K);
m = size(K,1);

non_zeros = diag(SK) >= 1e-11;
[C, ~, ~] = unique(diag(SK));
if ((sum(non_zeros)>=(m-1)) && (length(C)==m))
    % Only VX has degrees of freedom
    VXfixed = VK(:,1:m);
    VX0 = VK(:,(m+1):end);

    opts.record = 0;
    opts.mxitr  = 10000;
    opts.xtol = 1e-5;
    opts.gtol = 1e-5;
    opts.ftol = 1e-8;
    opts.tau = 1e-3;
    reg = 20e2;

    [VXout, ~] = OptStiefelGBB(VX0, @(x)l1_norm_ind(x, VXfixed, reg), opts);
    VKopt = VK;
    VKopt(:,(m+1):end) = VXout;
    
    disp('Original matrix : ');
    disp(UK*SK*VK');

    disp('Reconstructed : ');
    disp(UK*SK*VKopt');

    disp('Original V : ');
    disp(VK);

    disp('Optimized V : ');
    disp(VKopt);
else
    % UX and VX both have degrees of freedom
    disp('Implement!');
    return
end
    function [F,G] = l1_norm_ind(x, x_fixed, reg)
        F = sum(abs(x), 'all');
        ind = x_fixed' * x; % p x m, should be zero
        F = F + reg * sum(ind.^2, 'all'); % Add constraint violation with lagrange multiplier

        G = sign(x) + 2*reg*x_fixed*ind;
    end
end