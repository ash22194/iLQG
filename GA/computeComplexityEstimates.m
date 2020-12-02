function err_compute = computeComplexityEstimates(sys, p, s)
%% 
% p is m x 2 matrix where pi1 denotes the parent to input i
%                         pi2 denotes the child ID for input i
% If two inputs have the same non-zero parent and the same child ID, then they are coupled
% s is m x n matrix where sij = 1 if input i is dependent on state j

%% 

[c, c_eq] = constraints(p, s);
if (any(c_eq~=0) || any(c > 0))
    err_compute = 1;
else
    p = round(p);
    s = logical(round(s));
    assert(any(p(:,1)==0), 'Root of the tree must be 0');
    assert(all(p(:,1)~=linspace(1,sys.U_DIMS,sys.U_DIMS)'), 'No input can be its own parent');

    p_ = [linspace(1, sys.U_DIMS, sys.U_DIMS)', p];
    % Build tree
    action_tree = {};
    queue = {{0; -1}};

    while(~isempty(queue))
        curr_node = queue{1}; % pop the top element
        curr_parent = curr_node{2};
        curr_node = curr_node{1};
        queue = queue(2:end);

        children = p_(any(p_(:,2)==curr_node', 2), [1, 3]); % Find inputs that are children to curr_parent
        childID = unique(children(:,2));                    % Find unique childIDs, inputs are coupled if child IDs are same
        curr_children = cell(length(childID), 1);
        for ii=1:1:length(childID)
            curr_children{ii} = children(children(:,2)==childID(ii), 1);
            assert(all(s(curr_children{ii}, :) == prod(s(curr_children{ii}, :), 1), 'all'), ...
                   'Coupled inputs must have same state assignment');
            queue{end+1} = {curr_children{ii}; curr_node};
        end

        if (curr_node~=0)
            curr_state = s(curr_node(1), :);
        else
            curr_state = zeros(1, sys.X_DIMS);
        end
        action_tree(end+1, :) = {curr_node, curr_state, ...
                                 curr_parent, curr_children};
    end
    
%     err_compute = -sys.U_DIMS*prod(sys.numPoints);
    err_compute = 0;
    while (~isempty(action_tree))
        % Find leaf nodes
        leaf_node_ids = cellfun(@(x) isempty(x), action_tree(:,end));
        leaf_nodes = action_tree(leaf_node_ids, :);
%         leaf_node_ids = find(leaf_node_ids);

        if (size(leaf_nodes, 1)==1 && all(leaf_nodes{1,1}==0))
            break;
        end
        
        for ii=1:1:size(leaf_nodes,1)
            NS = sys.num_points(logical(leaf_nodes{ii,2}));
            NA = sys.num_action_samples(leaf_nodes{ii,1});
            M = sys.max_iter;
            MP = sys.max_policy_iter;
            interp_complexity = 2^length(NS); % Gridded Interpolation
            step_complexity = 4 * length(NS); % RK4 integration
            sample_complexity = 1; % randomly sample action
            action_update_complexity = 2;
            subpolicy_eval_compute = prod(NS) * (MP * interp_complexity + step_complexity);
            subpolicy_update_compute = prod(NS) * prod(NA) ...
                                       * (sample_complexity + step_complexity ...
                                          + interp_complexity + action_update_complexity);
                                      
            err_compute = err_compute + M * (subpolicy_eval_compute + subpolicy_update_compute);
            
            parent_input = leaf_nodes{ii, end-1};
            parent_node_id = find(cellfun(@(x) isempty(setdiff(x, parent_input)) && isempty(setdiff(parent_input, x)),...
                                  action_tree(:, 1)));
            assert(length(parent_node_id)==1, 'Invalid Parent Node ID');
            parent_node = action_tree(parent_node_id, :);

            assert(all(~(parent_node{2} .* leaf_nodes{ii, 2})), 'Invalid state overlap');
            parent_node{2} = parent_node{2} + leaf_nodes{ii, 2};

            % Find and delete the leaf_node in the list of children
            children_list = parent_node{end};
            childID = cellfun(@(x) isempty(setdiff(x, leaf_nodes{ii, 1})) ...
                                   && isempty(setdiff(leaf_nodes{ii, 1}, x)), ...
                            children_list);
            children_list(childID) = [];
            parent_node{end} = children_list;

            action_tree(parent_node_id, :) = parent_node;
        end
        action_tree(leaf_node_ids, :) = [];
    end
    NS = sys.num_points;
    NA = sys.num_action_samples;
    M = sys.max_iter;
    MP = sys.max_policy_iter;
    interp_complexity = 2^length(NS); % Gridded Interpolation
    step_complexity = 4 * length(NS); % RK4 integration
    sample_complexity = 1; % randomly sample action
    action_update_complexity = 2;
    policy_eval_compute = prod(NS) * (MP * interp_complexity + step_complexity);
    policy_update_compute = prod(NS) * prod(NA) ...
                            * (sample_complexity + step_complexity ...
                               + interp_complexity + action_update_complexity);
    joint_compute = M * (policy_eval_compute + policy_update_compute);
    
    err_compute = err_compute / joint_compute;
end

end