function encoding = encode_tree(action_tree)
    
    X_DIMS = sum(cellfun(@(x) sum(x), action_tree(:,2)));
    U_DIMS = sum(cellfun(@(x) length(x), action_tree(:,1))) - 1;
    state_assignment = zeros(U_DIMS, X_DIMS);
    % action tree = {input_tuple, state_dims, }
    % Conver the children to single values
    for ii=1:1:size(action_tree, 1)
        action_tree{ii, end} = cellfun(@(x) vec_to_scalar(x, U_DIMS + 1), action_tree{ii, end});
        if (all(action_tree{ii, 1}~=0))
            state_assignment(action_tree{ii, 1}, :) = ones(length(action_tree{ii, 1}), 1) ...
                                                        * double(action_tree{ii, 2});
        end
    end
    % Convert all input tuples to single values
    action_tree(:,1) = cellfun(@(x) vec_to_scalar(x, U_DIMS + 1), action_tree(:,1), ...
                               'UniformOutput', false);
    % Convert the parents to single values
    action_tree(:,3) = cellfun(@(x) vec_to_scalar(x, U_DIMS + 1), action_tree(:,3), ...
                               'UniformOutput', false);
        
    % Convert action tree to prufer code
    prufer_seq = [];
    while (size(action_tree, 1) > 2)
        leaf_node_ids = cellfun(@(x) isempty(x), action_tree(:,end));
        leaf_nodes = cell2mat(action_tree(leaf_node_ids, 1));
        leaf_node_parents = cell2mat(action_tree(leaf_node_ids, 3));
        leaf_node_ids = find(leaf_node_ids);

        [smallest_leaf_node, smallest_leaf_node_id] = min(leaf_nodes);
        smallest_leaf_node_parent = leaf_node_parents(smallest_leaf_node_id);
        smallest_leaf_node_parent_id = find(cellfun(@(x) smallest_leaf_node_parent==x, action_tree(:,1)));

        prufer_seq = [prufer_seq, smallest_leaf_node_parent];
        children = action_tree{smallest_leaf_node_parent_id, end};
        children_to_remove = children==smallest_leaf_node;
        children(children_to_remove) = [];
        action_tree{smallest_leaf_node_parent_id, end} = children;
        action_tree(leaf_node_ids(smallest_leaf_node_id),:) = [];
    end
    prufer_seq = [prufer_seq, ...
                  -1*ones(1, U_DIMS - length(prufer_seq) - 1)];
    assert(length(prufer_seq) == (U_DIMS - 1), 'Check sequence length');
    
    encoding = 48 + [prufer_seq, reshape(state_assignment, 1, U_DIMS*X_DIMS)];
end

function scalar = vec_to_scalar(vec, modulo)
    vec = vec(:);
    scalar = sum(sort(vec) .* (modulo.^((0:1:(size(vec,1)-1))')));
end