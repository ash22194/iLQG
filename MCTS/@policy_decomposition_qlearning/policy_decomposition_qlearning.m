classdef policy_decomposition_qlearning < handle
    properties
        sys
        p
        s
        action_tree
        decomposition_key
        decomposition_key_numeric
        parent
        lqr_measure
        compute_fraction
        measure
        visit_count
        epsilon
        nodelist
        childnodes
        childnodes_encoding
        childnodes_unexpanded
        childnodes_subtree_unexplored
        childnodes_measure_prediction
        childnodes_measure_prediction_version
        childnode_maxcount
        childnode_best
        childnode_worst
        subtree_unexplored
    end
    
    methods 
        
        function de = policy_decomposition_qlearning(varargin)
            %% Description
            % Each node stores a list of child nodes that are objects of
            % type policy_decomposition_qlearning when expanded and string
            % encodings when unexpanded.
            % Only the pseudo_inputs that are leaves are split to create
            % new decompositions
            % Copies of expanded nodes are stored in a dictionary with string
            % encodings as keys
            % A neural network is trained to add progressive bias to the
            % selection strategy
            
            % Constructor
            sys = varargin{1};
            if (~isfield(sys, 'fxfu_func'))
                sys.fxfu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
            end
            
            if (~isfield(sys, 'err_lqr_func'))
                fxfu = sys.fxfu_func(sys.l_point, sys.u0);
                A = fxfu(:, 1:sys.X_DIMS);
                B = fxfu(:, (1+sys.X_DIMS):end);
                lambda_ = (1 - sys.gamma_) / sys.dt;
                [~, S_j, ~] = lqr(A - (lambda_ / 2) * eye(size(A,1)), B, sys.Q, sys.R, zeros(size(B)));
                
                sys.A = A;
                sys.B = B;
                sys.lambda_ = lambda_;
                sys.S = sym('S', [sys.X_DIMS, sys.X_DIMS]);
                sys.S = tril(sys.S,0) + tril(sys.S,-1).';
                sys.a = sym('a', [sys.X_DIMS, 1]);
                sys.b = sym('b', [sys.X_DIMS, 1]);
                err_lqr = sys.x.'*sys.S*sys.x - sys.x.'*S_j*sys.x;
                for ii=1:1:sys.X_DIMS
                    err_lqr = int(err_lqr, sys.x(ii), [sys.a(ii), sys.b(ii)]);
                end
                sys.err_lqr_func = matlabFunction(err_lqr, 'Vars', {sys.S, sys.a, sys.b});
                sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1), 1);
            end
            
            if (~isfield(sys, 'measure_func'))
                sys.measure_func = @(err_lqr, err_compute) err_lqr;
                % sys.measure_func = @(err_lqr, err_compute) err_lqr * err_compute;
                % sys.measure_func = @(err_lqr, err_compute) (1 - exp(-err_lqr)) * err_compute;
            end
            
            if (nargin==3)
                action_tree = varargin{2};
                parent      = varargin{3};
                
                % Construct p and s
                X_DIMS = sum(cellfun(@(x) sum(x), action_tree(:,2)));
                U_DIMS = sum(cellfun(@(x) length(x), action_tree(:,1))) - 1;
                assert(X_DIMS==sys.X_DIMS, sprintf('Check X_DIMS: (%d), (%d)', X_DIMS, sys.X_DIMS));
                assert(U_DIMS==sys.U_DIMS, sprintf('Check U_DIMS: (%d), (%d)', U_DIMS, sys.U_DIMS));
                
                p = zeros(U_DIMS, 2);
                s = zeros(U_DIMS, X_DIMS);
                
                for jj=1:1:size(action_tree, 1)
                    curr_node = action_tree{jj, 1};
                    if (all(curr_node==0))
                        continue;
                    end
                    curr_state = action_tree{jj, 2};
                    curr_parent = action_tree{jj, 3};
                    
                    parent_node = cellfun(@(x) fastsetequal(x, curr_parent), action_tree(:,1));
                    assert(sum(parent_node)==1, 'Incorrect parent node for child!');
                    
                    parent_node = action_tree(parent_node, :);
                    children_list = parent_node{1, end};
                    child_id = find(cellfun(@(x) fastsetequal(x, curr_node), children_list), 1);
                    assert(~isempty(child_id), 'Child not found!');
                    
                    p(curr_node, 1) = curr_parent(1);
                    p(curr_node, 2) = child_id;
                    
                    s(curr_node, curr_state) = 1;
                end
                
                assert(any(p(:,1)==0), 'Root of the tree must be 0');
                assert(all(p(:,1)~=linspace(1, U_DIMS, U_DIMS)'), 'No input can be its own parent');
                
            elseif (nargin==4)
                p      = varargin{2};
                s      = varargin{3};
                parent = varargin{4};
                assert(size(p,1) == size(s,1), 'Not a valid decomposition');
                
                % Build action tree
                X_DIMS = size(s, 2);
                U_DIMS = size(s, 1);
                assert(X_DIMS==sys.X_DIMS, sprintf('Check X_DIMS: (%d), (%d)', X_DIMS, sys.X_DIMS));
                assert(U_DIMS==sys.U_DIMS, sprintf('Check U_DIMS: (%d), (%d)', U_DIMS, sys.U_DIMS));
                
                p = round(p);
                s = logical(round(s));
                assert(any(p(:,1)==0), 'Root of the tree must be 0');
                assert(all(p(:,1)~=linspace(1, U_DIMS, U_DIMS)'), 'No input can be its own parent');

                p_ = [linspace(1, U_DIMS, U_DIMS)', p];
                action_tree = {};
                queue = {{0; -1}};

                while(~isempty(queue))
                    curr_node = queue{1}; % pop the top element
                    curr_parent = curr_node{2};
                    curr_node = curr_node{1};
                    if (curr_node~=0)
                        curr_state = s(curr_node(1), :);
                    else
                        curr_state = zeros(1, X_DIMS);
                    end

                    queue = queue(2:end);
                    children = p_(any(p_(:,2)==curr_node', 2), [1, 3]); % Find inputs that are children to curr_parent
                    childID = unique(children(:,2));                    % Find unique childIDs, inputs are coupled if child IDs are same
                    curr_children = cell(length(childID), 1);
                    for jj=1:1:length(childID)
                        curr_children{jj} = children(children(:,2)==childID(jj), 1);
                        assert(all(s(curr_children{jj}, :) == prod(s(curr_children{jj}, :), 1), 'all'), ...
                               'Coupled inputs must have same state assignment');

                        queue{end+1} = {curr_children{jj}; curr_node};
                    end
                    action_tree(end+1, :) = {curr_node, curr_state, curr_parent, curr_children};
                end
            else
                disp('Number of inputs must be 3 or 4');
                return;
            end
            
            de.sys = sys;
            de.p = p;
            de.s = s;
            de.action_tree = action_tree;
            encoding = encode_binary_adj(de.action_tree);
            de.decomposition_key_numeric = encoding;
%             [A_, B_] = pre_process_dyn(sys, de.action_tree);
%             de.decomposition_key_numeric = [encoding, ...
%                                             reshape(A_, 1, sys.X_DIMS^2), ...
%                                             reshape(B_, 1, sys.U_DIMS * sys.X_DIMS)];
            
            de.decomposition_key = fastint2str(encoding);
            de.parent = parent;
            if (isa(parent, 'policy_decomposition_qlearning'))
                de.nodelist = parent.nodelist;
            else
                de.nodelist = containers.Map();
            end
            assert(~de.nodelist.isKey(de.decomposition_key), 'Node already exists!!');
            de.compute_measure();
            de.childnodes = [];
            de.childnodes_encoding = [];
            de.childnodes_unexpanded = [];
            de.childnodes_subtree_unexplored = [];
            de.childnodes_measure_prediction = [];
            de.childnodes_measure_prediction_version = -1;
            de.childnode_best = 0;
            de.childnode_worst = 0;
            de.visit_count = 0;
            de.epsilon = sys.epsilon;
            de.subtree_unexplored = true;
            
            % Calculate maximum children you can get from a node
            % Children can be created by splitting pseudo-inputs
            childnode_maxcount = 0;
            for jj=1:1:size(action_tree, 1)
                num_inputs = length(action_tree{jj, 1});
                num_states = sum(action_tree{jj, 2});
                is_leaf = isempty(action_tree{jj, end});
                
                num_input2groups = distr_n_into_r(num_inputs, 2) / 2; % Number of non-empty 2-groups
                if (is_leaf) % Is a leaf?
                    % Add a cascade
                    num_state2groups = distr_n_into_r_with_k_nonempty(num_states, 2, 1);
                    % Decouple the pseudo-inputs
                    num_state2groups_nonempty = distr_n_into_r(num_states, 2);
                    
                    decomp_count = num_input2groups * (2 * num_state2groups + num_state2groups_nonempty);
                    childnode_maxcount = childnode_maxcount + decomp_count;
                end
            end
            de.childnode_maxcount = childnode_maxcount;
            
            % Add node to nodelist
            de.nodelist(de.decomposition_key) = de;
        end
        
        function par = get_parent(de)
            par = de.parent;
        end
        
        function child = expand_best_child(de, determ, predictor)
            
            if isempty(de.childnodes)
                de.create_possible_children(predictor);
                if (isempty(de.childnodes))
                    de.visit_count = de.visit_count + 1;
                    de.subtree_unexplored = false;
                    child = [];
                    return;
                end
            end
            
            de.visit_count = de.visit_count + 1;
            
            % Update predictions if the predictor has been re-trained
            if (~(de.childnodes_measure_prediction_version == predictor.version))
                if (de.childnode_maxcount > 0)
                    de.childnodes_measure_prediction = predictor.predict(de.childnodes_encoding);
                end
                de.childnodes_measure_prediction_version = predictor.version;
            end

            if (sum(de.childnodes_subtree_unexplored)==0)
                de.subtree_unexplored = false;
                de.visit_count = de.visit_count - 1;
                child = de.parent;
                return;
            end

            % Use epsilon-greedy strategy to navigate
            explore = false;
            if (rand() <= de.epsilon)
                % random exploration
                explore = true;
            end

            is_valid = false;
            while(~is_valid)
                if (explore)
                    child_id = randi([1, sum(de.childnodes_subtree_unexplored)], 1);
                else
                    childnodes_measure = de.childnodes_measure_prediction;
                    childnodes_measure(~de.childnodes_unexpanded) ...
                        = cellfun(@(x) x.measure, de.childnodes(~de.childnodes_unexpanded));
                    [child_measure_min, child_id] = min(childnodes_measure);
                    if (~determ)
                        best_children = find(childnodes_measure==child_measure_min);
                        child_id = best_children(randi([1, length(best_children)], 1));
                    end
                end
                child_id = find(de.childnodes_subtree_unexplored, child_id);
                child_id = child_id(end);
                child = de.childnodes{child_id};
                if (~isa(child, 'policy_decomposition_qlearning'))
                    child_encoding = fastint2str(de.childnodes_encoding(1,:,1,child_id));
                    if (de.nodelist.isKey(child_encoding))
                        child = de.nodelist(child_encoding);
                        child.parent = de;
                    else
                        child = policy_decomposition_qlearning(de.sys, ...
                                                               child, ...
                                                               de);                   
                    end
                    de.childnodes{child_id} = child;
                    assert(de.childnodes_unexpanded(child_id));
                    de.childnodes_unexpanded(child_id) = false;
                end
                if (~child.subtree_unexplored)
                    de.childnodes_subtree_unexplored(child_id) = false;
                else
                    is_valid = true;
                end
                if (sum(de.childnodes_subtree_unexplored)==0)
                    de.subtree_unexplored = false;
                    [~, child_id] = min(cellfun(@(x) x.measure, de.childnodes)); % use the best found values
                    child = de.childnodes{child_id};
                    break;
                end
            end

            min_measure = child.measure;
            child_ = child;
            while (true)
                if (min_measure < de.measure)
                    de.measure = min_measure;
                    de.childnode_best = child_;
                    child_ = de;
                    de = de.parent;
                else
                    if ((de.childnode_best==0) || (min_measure < de.childnode_best.measure))
                        de.childnode_best = child_;
                    end
                    break;
                end
            end
        end
        
        function visualize_subtree(de)
            
            % Current node is the root
            nodes = [0];
            lqr_measures = [de.lqr_measure];
            compute_fractions = [de.compute_fraction];
            queue = {de};
            parent_count = 1;
            
            while(~isempty(queue))
                de_ = queue{1};
                queue = queue(2:end,:);
                if (isempty(de_.childnodes))
                    continue;
                end
                nodes = cat(1, nodes, parent_count*ones(size(de_.childnodes, 1), 1));
                lqr_measures = cat(1, lqr_measures, cellfun(@(x) x.lqr_measure, de_.childnodes));
                compute_fractions = cat(1, compute_fractions, cellfun(@(x) x.compute_fraction, de_.childnodes));
                
                queue = cat(1, queue, de_.childnodes);
                parent_count = parent_count + 1;
            end
            
            [x, y] = treelayout(nodes');
            treeplot(nodes'); hold on;
            text(x, y+0.01, num2str(lqr_measures, '%.3f'), 'FontSize', 5);
            text(x, y-0.01, num2str(compute_fractions, '%.3f'), 'FontSize', 5);
            hold off;
        end
        
        function [num, num_unique] = count_subtree_nodes(de)
            
            num = 0;
            num_unique = 0;
            unique_nodes = containers.Map();
            queue = {};
            if (de.visit_count > 0)
                queue = cat(1, queue, {de});
            end
            while (~isempty(queue))
                first_node = queue{1};
                queue = queue(2:end, :);
                if (~unique_nodes.isKey(first_node.decomposition_key))
                    unique_nodes(first_node.decomposition_key) = 1;
                    num_unique = num_unique + 1;
                end
                num = num + 1;
                
                if (~isempty(de.childnodes))
                    visited_children = cellfun(@(x) x.visit_count > 0, de.childnodes);
                    queue = cat(1, queue, de.childnodes(visited_children));
                end
            end
        end
    end
end