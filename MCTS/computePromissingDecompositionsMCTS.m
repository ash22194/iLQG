clear;
close all;
clc;

%% 

restoredefaultpath();
system_name = 'manipulator3dof';
addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
addpath('utils');
load(strcat('data/', system_name, 'System.mat'), 'sys');
sys.measure_func = @(err_lqr, err_compute) (1 - exp(-err_lqr)) .* err_compute;

max_mcts_iter = inf;
max_mcts_time = 100;
deterministic = false;

%% Alg 1)

mctree1 = policy_decomposition(sys, ...
                               [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)], ...
                               ones(sys.U_DIMS, sys.X_DIMS), ...
                               -1);
tic;
run_mcts(mctree1, max_mcts_iter, max_mcts_time, deterministic);
time_mcts1 = toc;

%% Alg 2)

mctree2 = policy_decomposition_wrep(sys, ...
                                   [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)], ...
                                   ones(sys.U_DIMS, sys.X_DIMS), ...
                                   -1);
tic;
run_mcts(mctree2, max_mcts_iter, max_mcts_time, deterministic);
time_mcts2 = toc;

%% Alg 3)

mctree3 = policy_decomposition_worep(sys, ...
                                   [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)], ...
                                   ones(sys.U_DIMS, sys.X_DIMS), ...
                                   {});
tic;
run_mcts(mctree3, max_mcts_iter, max_mcts_time, deterministic);
time_mcts3 = toc;

%% Extract Promissing Decompositions
num_to_extract = 10;

% Alg 1
root = mctree1;

children = root.childnodes(cellfun(@(x) isa(root.nodelist(x), class(root)), root.childnodes));
children = cellfun(@(x) root.nodelist(x), children, 'UniformOutput', false);
children_measure = cellfun(@(x) x.measure, children);
[~, children_measure_order] = sort(children_measure);
children = children(children_measure_order, :);

num_extracted = 0;
best_children1 = cell(0, 1);
while (num_extracted < num_to_extract && ~isempty(children))
    best_children1 = cat(1, best_children1, children(1,:));
    while (abs(best_children1{end, 1}.measure ...
                - sys.measure_func(best_children1{end, 1}.lqr_measure, ...
                                   best_children1{end, 1}.compute_fraction)) > eps)
         best_children1{end, 1} = best_children1{end, 1}.childnode_best;
    end
    assert(best_children1{end, 1}.measure == children{1, 1}.measure, 'Check the best child lookup');
    
    if (any(cellfun(@(x) strcmp(x.decomposition_key, ...
                                best_children1{end, 1}.decomposition_key), ...
                    best_children1(1:(end-1),:))))
        best_children1 = best_children1(1:(end-1), :);
    else
        num_extracted = num_extracted + 1;
    end
    children = children(2:end,:);
end

% Remove decompositions that are repeated
best_children_decomposition_id1 = cellfun(@(x) [reshape(x.p, 1, 2*sys.U_DIMS), reshape(x.s, 1, sys.U_DIMS*sys.X_DIMS)], ...
                                          best_children1, 'UniformOutput', false);
best_children_decomposition_id1 = cell2mat(best_children_decomposition_id1);

best_children_lqr_measure1 = cellfun(@(x) x.lqr_measure, best_children1);
best_children_compute_fraction1 = cellfun(@(x) x.compute_fraction, best_children1);

assert(all(abs(cellfun(@(x) x.measure, best_children1)...
               - sys.measure_func(best_children_lqr_measure1, best_children_compute_fraction1))...
            < eps), 'Check measure1 calculation!');

% Alg 2
root = mctree2;

children = root.childnodes(cellfun(@(x) isa(x, class(root)), root.childnodes));
children_measure = cellfun(@(x) x.measure, children);
[~, children_measure_order] = sort(children_measure);
children = children(children_measure_order, :);

num_extracted = 0;
best_children2 = cell(0, 1);
while (num_extracted < num_to_extract && ~isempty(children))
    best_children2 = cat(1, best_children2, children(1,:));
    while (abs(best_children2{end, 1}.measure ...
               - sys.measure_func(best_children2{end, 1}.lqr_measure, ...
                                  best_children2{end, 1}.compute_fraction)) > eps)
         best_children2{end, 1} = best_children2{end, 1}.childnode_best;
    end
    assert(best_children2{end, 1}.measure == children{1, 1}.measure, 'Check the best child lookup');
    
    if (any(cellfun(@(x) strcmp(x.decomposition_key, ...
                                best_children2{end, 1}.decomposition_key), ...
                    best_children2(1:(end-1),:))))
        best_children2 = best_children2(1:(end-1), :);
    else
        num_extracted = num_extracted + 1;
    end
    children = children(2:end,:);
end

% Remove decompositions that are repeated
best_children_decomposition_id2 = cellfun(@(x) [reshape(x.p, 1, 2*sys.U_DIMS), reshape(x.s, 1, sys.U_DIMS*sys.X_DIMS)], ...
                                          best_children2, 'UniformOutput', false);
best_children_decomposition_id2 = cell2mat(best_children_decomposition_id2);

best_children_lqr_measure2 = cellfun(@(x) x.lqr_measure, best_children2);
best_children_compute_fraction2 = cellfun(@(x) x.compute_fraction, best_children2);

assert(all(abs(cellfun(@(x) x.measure, best_children2)...
               - sys.measure_func(best_children_lqr_measure2, best_children_compute_fraction2))...
            < eps), 'Check measure2 calculation!');

% Alg 3
root = mctree3;

children = root.childnodes(cellfun(@(x) isa(x, class(root)), root.childnodes));
children_measure = cellfun(@(x) x.measure, children);
[~, children_measure_order] = sort(children_measure);
children = children(children_measure_order, :);

num_extracted = 0;
best_children3 = cell(0, 1);
while (num_extracted < num_to_extract && ~isempty(children))
    best_children3 = cat(1, best_children3, children(1,:));
    while (abs(best_children3{end, 1}.measure ...
               - sys.measure_func(best_children3{end, 1}.lqr_measure, ...
                                  best_children3{end, 1}.compute_fraction)) > eps)
         best_children3{end, 1} = best_children3{end, 1}.childnode_best;
    end
    assert(best_children3{end, 1}.measure == children{1, 1}.measure, 'Check the best child lookup');
    
    if (any(cellfun(@(x) strcmp(x.decomposition_key, ...
                                best_children3{end, 1}.decomposition_key), ...
                    best_children3(1:(end-1),:))))
        best_children3 = best_children3(1:(end-1), :);
    else
        num_extracted = num_extracted + 1;
    end
    children = children(2:end,:);
end

% Remove decompositions that are repeated
best_children_decomposition_id3 = cellfun(@(x) [reshape(x.p, 1, 2*sys.U_DIMS), reshape(x.s, 1, sys.U_DIMS*sys.X_DIMS)], ...
                                          best_children3, 'UniformOutput', false);
best_children_decomposition_id3 = cell2mat(best_children_decomposition_id3);

best_children_lqr_measure3 = cellfun(@(x) x.lqr_measure, best_children3);
best_children_compute_fraction3 = cellfun(@(x) x.compute_fraction, best_children3);

assert(all(abs(cellfun(@(x) x.measure, best_children3)...
               - sys.measure_func(best_children_lqr_measure3, best_children_compute_fraction3))...
            < eps), 'Check measure3 calculation!');

save(strcat('data/', system_name, '_MCTS_', num2str(max_mcts_time), '_run5.mat'), 'sys', ...
            'best_children_decomposition_id1', 'best_children_lqr_measure1', 'best_children_compute_fraction1', ...
            'best_children_decomposition_id2', 'best_children_lqr_measure2', 'best_children_compute_fraction2', ...
            'best_children_decomposition_id3', 'best_children_lqr_measure3', 'best_children_compute_fraction3');
