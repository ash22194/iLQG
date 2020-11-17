function root_node = run_mcts(sys, max_iter, max_time, determ)
    
    root_node = policy_decomposition(sys, ...
                                     [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)], ...
                                     ones(sys.U_DIMS, sys.X_DIMS), ...
                                     -1);
    root_node.compute_measure();
    
    % Expand nodes till you hit a tree
    ii = 0;
    tic;
    while ((ii < max_iter) && (toc < max_time) && (~root_node.subtree_explored))
        ii = ii + 1;
        curr_node = root_node;
        while(isa(curr_node, 'policy_decomposition'))
            try
                curr_node = curr_node.expand_best_child(determ);
            catch ME
                disp('Something is wrong!');
            end
        end
        fprintf('Iter : %d\n', ii);
        fprintf('Best so far (%d):\np=', root_node.childnode_best.measure);
        disp(num2str(root_node.childnode_best.p));
        fprintf('s=');
        disp(num2str(root_node.childnode_best.s));
    end
end