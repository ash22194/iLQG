function [root_node, info] = run_mcts_nn(root_node, max_iter, max_time, determ, nn_params)
    
    root_node.compute_measure();
    data_dim = length(root_node.decomposition_key_numeric);
    
    if (isfield(nn_params, 'net'))
        net = nn_params.net;
        predictor.predict = @(x) net.predict(x);
    else
        net = {};
        predictor.predict = @(x) 1e-6*ones(size(x,4), 1);
    end
    predictor.version = 0;
    train_data = zeros(1, data_dim, 1, 0);
    
    % Expand nodes till you hit a tree
    ii = 0;
    curr_best_measure = 1;
    tic;
    while ((ii < max_iter) && (toc < max_time) && (root_node.subtree_unexplored))
        ii = ii + 1;
        curr_node = root_node;
        while(isa(curr_node, class(root_node)))
            try
                curr_node = curr_node.expand_best_child(determ, predictor);
            catch ME
                disp('Something is wrong!');
            end
        end
        
        num_nodes_explored = length(root_node.nodelist.values); % Assuming only expanded nodes are added to the nodelist
        train_again = (num_nodes_explored - size(train_data, 4)) >= nn_params.train_every;

        if (train_again)
            % Train with the new data
            train_data = reshape(cell2mat(cellfun(@(x) x.decomposition_key_numeric, ...
                                                 root_node.nodelist.values', ...
                                                 'UniformOutput', false))', ...
                                 1, data_dim, 1, num_nodes_explored);
%             train_labels = cellfun(@(x) root_node.sys.measure_func(x.lqr_measure, x.compute_fraction), ...
%                                    root_node.nodelist.values');
            train_labels = cellfun(@(x) x.measure, ...
                                   root_node.nodelist.values');
            predictor_labels_before = predictor.predict(train_data);
            err_before = (train_labels - predictor_labels_before);
            if (mean(err_before)>1e-2)
                if (~isa(net, 'SeriesNetwork'))
                    [net, info] = trainNetwork(train_data, train_labels, nn_params.layers, nn_params.options);
                else
                    [net, info] = trainNetwork(train_data, train_labels, net.Layers, nn_params.options);
                end
                predictor_labels_after = net.predict(train_data);
                err_after = (train_labels - predictor_labels_after);
                retrain_count = 0;
                while(mean(err_after) > 1e-2)
                    [net, info] = trainNetwork(train_data, train_labels, net.Layers, nn_params.options);
                    predictor_labels_after = net.predict(train_data);
                    err_after = (train_labels - predictor_labels_after);
                    retrain_count = retrain_count + 1;
                    disp(strcat('Retraining : ', num2str(retrain_count)));
                end

                predictor.version = predictor.version + 1;
                predictor.predict = @(x) net.predict(x);

                disp(strcat('Number of nodes : ', num2str(num_nodes_explored), ...
                            ', Err before : ', num2str(mean((err_before).^2, 1)), ...
                                ', after : ', num2str(mean((err_after).^2, 1))));
            else
                disp('Error below threshold, no training required!');
            end
        end
        
        if (root_node.childnode_best.measure < curr_best_measure)
            time_to_find_best = toc;
            curr_best_measure = root_node.childnode_best.measure;
            fprintf('Iter : %d, Nodes : %d, Best so far : %d\n', ...
                    ii, num_nodes_explored, curr_best_measure);
        end
    end
    
	num_nodes_explored = sum(cellfun(@(x) isa(x, class(root_node)), root_node.nodelist.values));
    
    fprintf('Best : %d, Time for Best : %d s, Number of nodes : %d\n', curr_best_measure, time_to_find_best, num_nodes_explored);
    
    info.best_measure = curr_best_measure;
    info.time_to_find_best = time_to_find_best;
    info.num_nodes_explored = num_nodes_explored;
end