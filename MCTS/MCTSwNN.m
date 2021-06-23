clear;
close all;
% clc;

%% 

restoredefaultpath();
system_name = 'manipulator4dof';
addpath('../iLQG/GA');
addpath('../iLQG/iLQG_boxQP/iLQG utilities/decomposition_count');
addpath(strcat('../iLQG/iLQG_boxQP/new_systems/', system_name));
addpath('utils');
load(strcat('data/', system_name, 'System.mat'), 'sys');
sys.measure_func = @(err_lqr, err_compute) (1 - exp(-err_lqr)) .* err_compute;
sys.X_DIMS_MASKED = eye(sys.X_DIMS);
                 
max_mcts_iter = inf;
max_mcts_time = 800;
num_runs = 5;
deterministic = false;

%% Q Learning

% Define the neural network
hidden_layer_size = 128;
dropout_rate = 0.25;
input_size = sys.U_DIMS*(sys.U_DIMS + sys.X_DIMS);% + sys.X_DIMS * (sys.X_DIMS + sys.U_DIMS);
nn_params_ql.layers = [imageInputLayer([1, input_size, 1], 'Normalization','none', 'Name', 'in')
                    fullyConnectedLayer(hidden_layer_size, 'Name', 'fc1')
                    reluLayer('Name', 'act1')
                    dropoutLayer(dropout_rate, 'Name', 'drop1') % layer 1 
                    fullyConnectedLayer(hidden_layer_size, 'Name', 'fc2')
                    reluLayer('Name', 'act2')
                    dropoutLayer(dropout_rate, 'Name', 'drop2') % layer 2
%                     fullyConnectedLayer(hidden_layer_size, 'Name', 'fc3')
%                     reluLayer('Name', 'act3')
%                     dropoutLayer(dropout_rate, 'Name', 'drop3') % layer 3
                    fullyConnectedLayer(hidden_layer_size, 'Name', 'out_fc1')
                    reluLayer('Name', 'out_act1')
                    dropoutLayer(dropout_rate, 'Name', 'out_drop1') % layer 4
                    fullyConnectedLayer(1, 'Name', 'out')
                    sigmoidLayer('Name', 'out_act')
                    regressionLayer()];
nn_params_ql.options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MiniBatchSize', 128, ...
    'MaxEpochs',100, ...
    'Shuffle','every-epoch', ...
    'Verbose',false);
nn_params_ql.train_every = 1500;
sys.epsilon = 0.1;

% Pre-train the neural network
% num_pretrain_points = 100000;
% population = uniform_input_tree_sampling(sys, num_pretrain_points);
% [~, uid] = extract_decompositions_from_population(sys, population);
% population = population(uid, :);
% population_cell = num2cell(population, 2);
% population_fitness = cellfun(@(x) computeJointMeasure(sys, ...
%                                                       reshape(x(1:(2*sys.U_DIMS)), sys.U_DIMS, 2), ...
%                                                       reshape(x((1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS)), ...
%                              population_cell);
% train_data = cellfun(@(x) encode_binary_adj_ps(reshape(x(1:(2*sys.U_DIMS)), sys.U_DIMS, 2), ...
%                                                reshape(x((1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS)), ...
%                      population_cell, 'UniformOutput', false);
% train_data = reshape(cell2mat(train_data)',1,sys.U_DIMS*(sys.U_DIMS + sys.X_DIMS),1,length(uid));
% [nn_params.net, ~] = trainNetwork(train_data, population_fitness, nn_params.layers, nn_params.options);
% load(strcat('data/', system_name, 'Pretrained.mat'));
% nn_params_ql.net = net;

mctree_ql = cell(num_runs, 2);
time_mcts_ql = zeros(num_runs, 1);
for rr=1:1:num_runs
    mctree_ql{rr, 1} = policy_decomposition_qlearning(sys, ...
                                              [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)], ...
                                              ones(sys.U_DIMS, sys.X_DIMS), ...
                                              {});
    tic;
%     profile on;
    [~, mctree_ql{rr, 2}] = run_mcts_nn(mctree_ql{rr, 1}, max_mcts_iter, max_mcts_time, deterministic, nn_params_ql);
%     profile off;
%     profile viewer;
    time_mcts_ql(rr) = toc;
end

%% NN based heuristic modifications

% Define the neural network
% hidden_layer_size = 128;
% dropout_rate = 0.25;
% input_size = sys.U_DIMS*(sys.U_DIMS + sys.X_DIMS);% + sys.X_DIMS * (sys.X_DIMS + sys.U_DIMS);
% nn_params.layers = [imageInputLayer([1, input_size, 1], 'Normalization','none', 'Name', 'in')
%                     fullyConnectedLayer(hidden_layer_size, 'Name', 'fc1')
%                     reluLayer('Name', 'act1')
%                     dropoutLayer(dropout_rate, 'Name', 'drop1') % layer 1 
%                     fullyConnectedLayer(hidden_layer_size, 'Name', 'fc2')
%                     reluLayer('Name', 'act2')
%                     dropoutLayer(dropout_rate, 'Name', 'drop2') % layer 2
% %                     fullyConnectedLayer(hidden_layer_size, 'Name', 'fc3')
% %                     reluLayer('Name', 'act3')
% %                     dropoutLayer(dropout_rate, 'Name', 'drop3') % layer 3
%                     fullyConnectedLayer(hidden_layer_size, 'Name', 'out_fc1')
%                     reluLayer('Name', 'out_act1')
%                     dropoutLayer(dropout_rate, 'Name', 'out_drop1') % layer 4
%                     fullyConnectedLayer(1, 'Name', 'out')
%                     sigmoidLayer('Name', 'out_act')
%                     regressionLayer()];
% nn_params.options = trainingOptions('sgdm', ...
%     'InitialLearnRate',0.01, ...
%     'MiniBatchSize', 128, ...
%     'MaxEpochs',200, ...
%     'Shuffle','every-epoch', ...
%     'Verbose',false);
% nn_params.train_every = 1000;
% 
% % Pre-train the neural network
% % num_pretrain_points = 100000;
% % population = uniform_input_tree_sampling(sys, num_pretrain_points);
% % [~, uid] = extract_decompositions_from_population(sys, population);
% % population = population(uid, :);
% % population_cell = num2cell(population, 2);
% % population_fitness = cellfun(@(x) computeJointMeasure(sys, ...
% %                                                       reshape(x(1:(2*sys.U_DIMS)), sys.U_DIMS, 2), ...
% %                                                       reshape(x((1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS)), ...
% %                              population_cell);
% % train_data = cellfun(@(x) encode_binary_adj_ps(reshape(x(1:(2*sys.U_DIMS)), sys.U_DIMS, 2), ...
% %                                                reshape(x((1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS)), ...
% %                      population_cell, 'UniformOutput', false);
% % train_data = reshape(cell2mat(train_data)',1,sys.U_DIMS*(sys.U_DIMS + sys.X_DIMS),1,length(uid));
% % [nn_params.net, ~] = trainNetwork(train_data, population_fitness, nn_params.layers, nn_params.options);
% % load(strcat('data/', system_name, 'Pretrained.mat'));
% % nn_params.net = net;
% 
% % Test progressive biasing and progressive unpruning
% sys.unpruning_init_fraction = 0.5;
% sys.unpruning_fraction = 0.1;
% sys.unpruning_base_fraction = 0.5;
% sys.unpruning_rate = 1.01;
% 
% mctree_pb = cell(num_runs, 2);
% time_mcts_pb = zeros(num_runs, 1);
% for rr=1:1:num_runs
%     mctree_pb{rr, 1} = policy_decomposition_pbpu(sys, ...
%                                               [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)], ...
%                                               ones(sys.U_DIMS, sys.X_DIMS), ...
%                                               {});
%     tic;
% %     profile on;
%     [~, mctree_pb{rr, 2}] = run_mcts_nn(mctree_pb{rr, 1}, max_mcts_iter, max_mcts_time, deterministic, nn_params);
% %     profile off;
% %     profile viewer;
%     time_mcts_pb(rr) = toc;
% end

%% Vanilla MCTS

mctree_vanila = cell(num_runs, 2);
time_mcts_vanila = zeros(num_runs, 1);
for rr=1:1:num_runs
    mctree_vanila{rr, 1} = policy_decomposition_leafexpand(sys, ...
                                              [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)], ...
                                              ones(sys.U_DIMS, sys.X_DIMS), ...
                                              {});
    tic;
    [~, mctree_vanila{rr, 2}] = run_mcts(mctree_vanila{rr, 1}, max_mcts_iter, max_mcts_time, deterministic);
    time_mcts_vanila(rr) = toc;
end

%% Extract Promissing Decompositions
num_to_extract = 10;

% Alg Q Learning
algql_stats = zeros(num_runs, 8);
algql_best_decompositions = cell(num_runs, 3);
for rr=1:1:num_runs
    root = mctree_ql{rr, 1};
    info = mctree_ql{rr, 2};

    children = root.childnodes(cellfun(@(x) isa(x, class(root)), root.childnodes));
    children_measure = cellfun(@(x) x.measure, children);
    [~, children_measure_order] = sort(children_measure);
    children = children(children_measure_order, :);

    num_extracted = 0;
    best_childrenql = cell(0, 1);
    while (num_extracted < num_to_extract && ~isempty(children))
        best_childrenql = cat(1, best_childrenql, children(1,:));
        while (abs(best_childrenql{end, 1}.measure ...
                   - sys.measure_func(best_childrenql{end, 1}.lqr_measure, ...
                                      best_childrenql{end, 1}.compute_fraction)) > eps)
             best_childrenql{end, 1} = best_childrenql{end, 1}.childnode_best;
        end
        assert(best_childrenql{end, 1}.measure <= children{1, 1}.measure, 'Check the best child lookup');

        if (any(cellfun(@(x) strcmp(x.decomposition_key, ...
                                    best_childrenql{end, 1}.decomposition_key), ...
                        best_childrenql(1:(end-1),:))))
            best_childrenql = best_childrenql(1:(end-1), :);
        else
            num_extracted = num_extracted + 1;
        end
        children = children(2:end,:);
    end

    % Remove decompositions that are repeated
    best_children_decomposition_idql = cellfun(@(x) [reshape(x.p, 1, 2*sys.U_DIMS), reshape(x.s, 1, sys.U_DIMS*sys.X_DIMS)], ...
                                              best_childrenql, 'UniformOutput', false);
    best_children_decomposition_idql = cell2mat(best_children_decomposition_idql);

    best_children_lqr_measureql = cellfun(@(x) x.lqr_measure, best_childrenql);
    best_children_compute_fractionql = cellfun(@(x) x.compute_fraction, best_childrenql);

    assert(all(abs(cellfun(@(x) x.measure, best_childrenql)...
                   - sys.measure_func(best_children_lqr_measureql, best_children_compute_fractionql))...
                < eps), 'Check measure ql calculation!');
    
    algql_stats(rr, 1) = info.best_measure;
    algql_stats(rr, 2) = mean(cellfun(@(x) x.measure, best_childrenql));
    algql_stats(rr, 3) = std(cellfun(@(x) x.measure, best_childrenql));
    
    algql_stats(rr, 4) = min(best_children_lqr_measureql);
    algql_stats(rr, 5) = mean(best_children_lqr_measureql);
    algql_stats(rr, 6) = std(best_children_lqr_measureql);
    
    algql_stats(rr, 7) = info.time_to_find_best;
    algql_stats(rr, 8) = info.num_nodes_explored;
    
    algql_best_decompositions{rr, 1} = best_children_decomposition_idql;
    algql_best_decompositions{rr, 2} = best_children_lqr_measureql;
    algql_best_decompositions{rr, 3} = best_children_compute_fractionql;
end

% % Alg Progressive Biasing
% algpb_stats = zeros(num_runs, 8);
% algpb_best_decompositions = cell(num_runs, 3);
% for rr=1:1:num_runs
%     root = mctree_pb{rr, 1};
%     info = mctree_pb{rr, 2};
% 
%     children = root.childnodes(cellfun(@(x) isa(x, class(root)), root.childnodes));
%     children_measure = cellfun(@(x) x.measure, children);
%     [~, children_measure_order] = sort(children_measure);
%     children = children(children_measure_order, :);
% 
%     num_extracted = 0;
%     best_childrenpb = cell(0, 1);
%     while (num_extracted < num_to_extract && ~isempty(children))
%         best_childrenpb = cat(1, best_childrenpb, children(1,:));
%         while (abs(best_childrenpb{end, 1}.measure ...
%                    - sys.measure_func(best_childrenpb{end, 1}.lqr_measure, ...
%                                       best_childrenpb{end, 1}.compute_fraction)) > eps)
%              best_childrenpb{end, 1} = best_childrenpb{end, 1}.childnode_best;
%         end
%         assert(best_childrenpb{end, 1}.measure <= children{1, 1}.measure, 'Check the best child lookup');
% 
%         if (any(cellfun(@(x) strcmp(x.decomposition_key, ...
%                                     best_childrenpb{end, 1}.decomposition_key), ...
%                         best_childrenpb(1:(end-1),:))))
%             best_childrenpb = best_childrenpb(1:(end-1), :);
%         else
%             num_extracted = num_extracted + 1;
%         end
%         children = children(2:end,:);
%     end
% 
%     % Remove decompositions that are repeated
%     best_children_decomposition_idpb = cellfun(@(x) [reshape(x.p, 1, 2*sys.U_DIMS), reshape(x.s, 1, sys.U_DIMS*sys.X_DIMS)], ...
%                                               best_childrenpb, 'UniformOutput', false);
%     best_children_decomposition_idpb = cell2mat(best_children_decomposition_idpb);
% 
%     best_children_lqr_measurepb = cellfun(@(x) x.lqr_measure, best_childrenpb);
%     best_children_compute_fractionpb = cellfun(@(x) x.compute_fraction, best_childrenpb);
% 
%     assert(all(abs(cellfun(@(x) x.measure, best_childrenpb)...
%                    - sys.measure_func(best_children_lqr_measurepb, best_children_compute_fractionpb))...
%                 < eps), 'Check measure pb calculation!');
%     
%     algpb_stats(rr, 1) = info.best_measure;
%     algpb_stats(rr, 2) = mean(cellfun(@(x) x.measure, best_childrenpb));
%     algpb_stats(rr, 3) = std(cellfun(@(x) x.measure, best_childrenpb));
%     
%     algpb_stats(rr, 4) = min(best_children_lqr_measurepb);
%     algpb_stats(rr, 5) = mean(best_children_lqr_measurepb);
%     algpb_stats(rr, 6) = std(best_children_lqr_measurepb);
%     
%     algpb_stats(rr, 7) = info.time_to_find_best;
%     algpb_stats(rr, 8) = info.num_nodes_explored;
%     
%     algpb_best_decompositions{rr, 1} = best_children_decomposition_idpb;
%     algpb_best_decompositions{rr, 2} = best_children_lqr_measurepb;
%     algpb_best_decompositions{rr, 3} = best_children_compute_fractionpb;
% end

% Alg vanila
algvanila_stats = zeros(num_runs, 8);
algvanila_best_decompositions = cell(num_runs, 3);
for rr=1:1:num_runs
    root = mctree_vanila{rr, 1};
    info = mctree_vanila{rr, 2};

    children = root.childnodes(cellfun(@(x) isa(x, class(root)), root.childnodes));
    children_measure = cellfun(@(x) x.measure, children);
    [~, children_measure_order] = sort(children_measure);
    children = children(children_measure_order, :);

    num_extracted = 0;
    best_childrenvanila = cell(0, 1);
    while (num_extracted < num_to_extract && ~isempty(children))
        best_childrenvanila = cat(1, best_childrenvanila, children(1,:));
        while (abs(best_childrenvanila{end, 1}.measure ...
                   - sys.measure_func(best_childrenvanila{end, 1}.lqr_measure, ...
                                      best_childrenvanila{end, 1}.compute_fraction)) > eps)
             best_childrenvanila{end, 1} = best_childrenvanila{end, 1}.childnode_best;
        end
        assert(best_childrenvanila{end, 1}.measure <= children{1, 1}.measure, 'Check the best child lookup');

        if (any(cellfun(@(x) strcmp(x.decomposition_key, ...
                                    best_childrenvanila{end, 1}.decomposition_key), ...
                        best_childrenvanila(1:(end-1),:))))
            best_childrenvanila = best_childrenvanila(1:(end-1), :);
        else
            num_extracted = num_extracted + 1;
        end
        children = children(2:end,:);
    end

    % Remove decompositions that are repeated
    best_children_decomposition_idvanila = cellfun(@(x) [reshape(x.p, 1, 2*sys.U_DIMS), reshape(x.s, 1, sys.U_DIMS*sys.X_DIMS)], ...
                                              best_childrenvanila, 'UniformOutput', false);
    best_children_decomposition_idvanila = cell2mat(best_children_decomposition_idvanila);

    best_children_lqr_measurevanila = cellfun(@(x) x.lqr_measure, best_childrenvanila);
    best_children_compute_fractionvanila = cellfun(@(x) x.compute_fraction, best_childrenvanila);

    assert(all(abs(cellfun(@(x) x.measure, best_childrenvanila)...
                   - sys.measure_func(best_children_lqr_measurevanila, best_children_compute_fractionvanila))...
                < eps), 'Check measure vanila calculation!');
    
    algvanila_stats(rr, 1) = info.best_measure;
    algvanila_stats(rr, 2) = mean(cellfun(@(x) x.measure, best_childrenvanila));
    algvanila_stats(rr, 3) = std(cellfun(@(x) x.measure, best_childrenvanila));
    
    algvanila_stats(rr, 4) = min(best_children_lqr_measurevanila);
    algvanila_stats(rr, 5) = mean(best_children_lqr_measurevanila);
    algvanila_stats(rr, 6) = std(best_children_lqr_measurevanila);
    
    algvanila_stats(rr, 7) = info.time_to_find_best;
    algvanila_stats(rr, 8) = info.num_nodes_explored;
    
    algvanila_best_decompositions{rr, 1} = best_children_decomposition_idvanila;
    algvanila_best_decompositions{rr, 2} = best_children_lqr_measurevanila;
    algvanila_best_decompositions{rr, 3} = best_children_compute_fractionvanila;
end

% save(strcat('data/final2/', system_name, '_MCTS_explqrobj', num2str(max_mcts_time), '.mat'), 'sys', ...
%      'algql_stats', 'algql_best_decompositions', ...
%      'algvanila_stats', 'algvanila_best_decompositions');