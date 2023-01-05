function [measure, err_lqr, err_compute] = computeJointMeasureTransformed(sys, p, s)
%% 
% p is m x 2 matrix where pi1 denotes the parent to input i
%                         pi2 denotes the child ID for input i
% If two inputs have the same non-zero parent and the same child ID, then they are coupled
% s is m x n matrix where sij = 1 if input i is dependent on state j

%% 
err_lqr_lowerbound = -1e-10;
[c, c_eq] = constraints(p, s);
if (any(c_eq~=0) || any(c > 0))
    err_lqr = inf;
    err_compute = 1;
    measure = sys.measure_func(err_lqr, err_compute);
    
else
    p = round(p);
    s = logical(round(s));
    
    if (isfield(sys, 'decompositionlist'))
        % Encode decomposition
        % Compute encoding for decomposition
        action_coupling = zeros(sys.U_DIMS, sys.U_DIMS);
        action_dependence = zeros(sys.U_DIMS, sys.U_DIMS);
        actions = 1:1:sys.U_DIMS;
        states = s;
        parents = p;
        while(~isempty(actions))
            acoupled_ = all((ones(length(actions), 1) * states(1,:)) ...
                            == states, 2) ...
                        & all((ones(length(actions), 1) * parents(1,:)) ...
                              == parents, 2);
            acoupled = actions(acoupled_);
            action_coupling(acoupled, acoupled) = 1;

            achildren = any(p(:,1) == acoupled, 2);
            action_dependence(acoupled, achildren) = 1;

            actions(acoupled_) = [];
            states(acoupled_, :) = [];
            parents(acoupled_, :) = [];
        end

        decomposition_key = [reshape(action_coupling, 1, sys.U_DIMS^2), ...
                             reshape(action_dependence, 1, sys.U_DIMS^2), ...
                             reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
        decomposition_key = fastint2str(decomposition_key);
        if (sys.decompositionlist.isKey(decomposition_key))
            measure = sys.decompositionlist(decomposition_key);
            err_lqr = measure(2);
            err_compute = measure(3);
            measure = measure(1);
            return;
        end
    end
    
    p_ = [linspace(1, sys.U_DIMS, sys.U_DIMS)', p];
    % Build tree
    action_tree = {};
    queue = {{0; -1}};
    action_coupling2 = zeros(sys.U_DIMS, sys.U_DIMS);
    action_dependence2 = zeros(sys.U_DIMS, sys.U_DIMS);
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

        u0_curr = -sys.Tv * sys.u0;
        K_curr = zeros(sys.U_DIMS, sys.X_DIMS);
        if (curr_node~=0)
            curr_state = s(curr_node(1), :);
            u0_curr(curr_node) = 0;
            action_coupling2(curr_node, curr_node) = 1;
            action_dependence2(curr_node, cell2mat(curr_children)) = 1;
        else
            curr_state = zeros(1, sys.X_DIMS);
        end
        action_tree(end+1, :) = {curr_node, curr_state, ...
                                 u0_curr, K_curr, ...
                                 curr_parent, {}, curr_children};
    end
    
    if (isfield(sys, 'decompositionlist'))
        assert(all(action_coupling == action_coupling2, 'all') && ...
               all(action_dependence == action_dependence2, 'all'), 'Check encoding scheme!');
    end

    err_compute = 0;
    while (~isempty(action_tree))
        % Find leaf nodes
        leaf_node_ids = cellfun(@(x) isempty(x), action_tree(:,end));
        leaf_nodes = action_tree(leaf_node_ids, :);

        if (size(leaf_nodes, 1)==1 && all(leaf_nodes{1,1}==0))
            K = leaf_nodes{1,4};
            try
                S = icare((sys.A - sys.B*K - sys.lambda_/2*eye(size(sys.A,1))), [], K'*sys.R*K + sys.Q, [], [], [], []);
            catch ME
                disp(ME.identifier);
                S = -eye(sys.X_DIMS);
            end
            
            eigS = eig(S);
            if (isempty(eigS) || any(~isreal(eigS)) || any(eigS < -1e-8))
                err_lqr = inf;
                err_compute = 1;
            else
                try
                    err_lqr = sys.err_lqr_func(S, sys.state_bounds(:,1)-sys.l_point, sys.state_bounds(:,2)-sys.l_point)/sys.da;
                    assert(err_lqr >= err_lqr_lowerbound, 'LQR error measure cannot be negative');
                catch
                    disp(ME.identifier);
                end

                NS = sys.num_points;
                NA = sys.num_action_samples;
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
                joint_compute = M * (subpolicy_eval_compute + subpolicy_update_compute);
                
                err_compute = err_compute / joint_compute;
            end
            err_lqr = abs(err_lqr);
            measure = sys.measure_func(err_lqr, err_compute);
            
            if ((measure == 0) && (err_lqr~=0))
                measure = err_lqr * err_compute;
            end
            
            if (isfield(sys, 'decompositionlist'))
                sys.decompositionlist(decomposition_key) = [measure, err_lqr, err_compute];
            end
            return;
        end
        
        for ii=1:1:size(leaf_nodes,1)
            u0_ = leaf_nodes{ii, 3};
            K = leaf_nodes{ii, 4};
            try
                fxfu = sys.fxfu_func(zeros(sys.X_DIMS, 1), u0_);
            catch
                disp('Error computing derivatives!');
            end
            
            A_ = fxfu(:,1:sys.X_DIMS);
            B_ = fxfu(:,(sys.X_DIMS+1):end);
            Q_ = sys.Q;
            R_ = sys.R;

            A_ = A_ - B_*K;
            A_ = A_(logical(leaf_nodes{ii, 2}), logical(leaf_nodes{ii, 2}));
            B_ = B_(logical(leaf_nodes{ii, 2}), leaf_nodes{ii, 1});
            Q_ = Q_ + K'*R_*K;
            Q_ = Q_(logical(leaf_nodes{ii, 2}), logical(leaf_nodes{ii, 2}));
            R_ = R_(leaf_nodes{ii, 1}, leaf_nodes{ii, 1});

            try
                [K_, ~, ~] = lqr(A_ - eye(size(A_, 1))*sys.lambda_/2, B_, ...
                                 Q_, R_, zeros(size(A_, 1), size(B_, 2)));
            catch ME
                err_lqr = inf;
                err_compute = 1;
                measure = sys.measure_func(err_lqr, err_compute);
                
                if (isfield(sys, 'decompositionlist'))
                    sys.decompositionlist(decomposition_key) = [measure, err_lqr, err_compute];
                end
                return;
            end
            K(leaf_nodes{ii, 1}, logical(leaf_nodes{ii, 2})) = K_;
            
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
            
            parent_input = leaf_nodes{ii, 5};
            parent_node_id = find(cellfun(@(x) isempty(setdiff(x, parent_input)) && isempty(setdiff(parent_input, x)),...
                                  action_tree(:, 1)));
            assert(length(parent_node_id)==1, 'Invalid Parent Node ID');
            parent_node = action_tree(parent_node_id, :);

            assert(all(~(parent_node{2} .* leaf_nodes{ii, 2})), 'Invalid state overlap');
            parent_node{2} = parent_node{2} + leaf_nodes{ii, 2};
            parent_node{3}(leaf_nodes{ii,1}) = leaf_nodes{ii, 3}(leaf_nodes{ii,1});
            for cch=1:1:size(leaf_nodes{ii,6}, 1)
                parent_node{3}(leaf_nodes{ii,6}{cch}) = leaf_nodes{ii, 3}(leaf_nodes{ii,6}{cch});
            end
            parent_node{4} = parent_node{4} + K;

            % Find and delete the leaf_node in the list of children
            children_list = parent_node{end};
            childID = cellfun(@(x) isempty(setdiff(x, leaf_nodes{ii, 1})) ...
                                   && isempty(setdiff(leaf_nodes{ii, 1}, x)), ...
                            children_list);
            children_list(childID) = [];
            parent_node{end} = children_list;
            parent_node{6} = cat(1, parent_node{6}, leaf_nodes{ii,1}, leaf_nodes{ii,6});

            action_tree(parent_node_id, :) = parent_node;
        end
        action_tree(leaf_node_ids, :) = [];
    end
end

end