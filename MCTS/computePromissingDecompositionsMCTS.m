clear;
close all;
clc;

%% 

restoredefaultpath();
system_name = 'manipulator3dof';
addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
addpath('utils');
load(strcat('data/', system_name, 'System.mat'), 'sys');
sys.measure_func = @(err_lqr, err_compute) (1 - exp(-err_lqr)) * err_compute;

num_mcts_runs = 1;
max_mcts_iter = inf;
max_mcts_time = 600;
deterministic = false;

mctree = cell(num_mcts_runs, 1);
tic;
for r=1:1:num_mcts_runs
    disp(strcat('Run : ', num2str(r)));
    mctree{r} = run_mcts(sys, max_mcts_iter, max_mcts_time, deterministic);
end
time_mcts = toc;

%% Extract Promissing Decompositions

num_to_extract = 10;
best_children = cell(num_to_extract, num_mcts_runs);
for r=1:1:num_mcts_runs
    
    root = mctree{r};
    child_measure = cellfun(@(x) x.measure, root.childnodes);
    [child_measure_sorted, child_measure_order] = sort(child_measure);
    best_children(:, r) = root.childnodes(child_measure_order(1:num_to_extract),:);
    
    for n=1:1:num_to_extract
        while (best_children{n, r}.measure ~= ...
                sys.measure_func(best_children{n, r}.lqr_measure, ...
                                 best_children{n, r}.compute_fraction))
             best_children{n, r} = best_children{n, r}.childnode_best;
        end
        assert(best_children{n, r}.measure == child_measure_sorted(n), 'Check the best child lookup');
    end
end

% Remove decompositions that are repeated
addpath('../GA');

best_children = reshape(best_children, num_mcts_runs * num_to_extract, 1);
best_children_decomposition_id = cellfun(@(x) [reshape(x.p, 1, 2*sys.U_DIMS), reshape(x.s, 1, sys.U_DIMS*sys.X_DIMS)], ...
                                          best_children, 'UniformOutput', false);
best_children_decomposition_id = cell2mat(best_children_decomposition_id);
[~, best_children_ids] = extract_decompositions_from_population(sys, best_children_decomposition_id);
best_children_decomposition_id = best_children_decomposition_id(best_children_ids,:);
best_children = best_children(best_children_ids,:);

[~, best_children_order] = sort(cellfun(@(x) x.measure, best_children));
best_children_decomposition_id = best_children_decomposition_id(1:min(num_to_extract, ...
                                                                      length(best_children)), :);
best_children = best_children(best_children_order(1:min(num_to_extract, ...
                                                        length(best_children))));
best_children_lqr_measure = cellfun(@(x) x.lqr_measure, best_children);
best_children_compute_fraction = cellfun(@(x) x.compute_fraction, best_children);

% save(strcat('data/', system_name, '_MCTS_', num2str(max_mcts_time), '.mat'), ...
%      'sys', 'best_children_decomposition_id', 'best_children_lqr_measure', 'best_children_compute_fraction');
