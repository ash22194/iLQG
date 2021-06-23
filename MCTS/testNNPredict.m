clear;
close all;
clc;

%%

% model = importKerasNetwork('data/manipulator4dof_elqr_predictor.h5');
% trainstat = load('data/manipulator4dof_trainstat.mat');
% load('data/manipulator4dofSystem.mat');
addpath('utils');

%% Test preprocess
u_x = [];
n = 3;
pseudo_inputs{1,1} = {{1, [2,3]}; {2, [1,3]}; {3, [1,2]}}; % r = 2
pseudo_inputs{2,1} = {{1, 2, 3}}; % r = 3

state_assignments{1} = {{1, [2,3]}; {2, [1,3]}; {3, [1,2]}; {[], [1,2,3]}};
state_assignments{2} = {{1, 2, 3}; {1, [2,3], []}; {2, [1,3], []}; {3, [1,2], []}; {[], [1,2,3], []}};

for r=1:1:size(pseudo_inputs, 1)
    for sub_r=1:1:size(pseudo_inputs{r}, 1)
        current_pseudo_inputs = pseudo_inputs{r}{sub_r};
        orders_of_inputs = perms(1:1:(r+1));
        for ori=1:1:size(orders_of_inputs,1)
            input_order = orders_of_inputs(ori,:);
            ordered_pseudo_input = current_pseudo_inputs(input_order);
            % Define action tree for purely cascaded
            p_casc = [zeros(n,1), ones(n,1)];
            for ori_=2:1:(r+1)
                p_casc(ordered_pseudo_input{ori_}, 1) = ordered_pseudo_input{ori_-1}(1);
            end

            for sub_s=1:1:size(state_assignments{r}, 1)
                current_state_assignment = state_assignments{r}{sub_s};
                orders_of_states = perms(1:1:(r+1));
                perms_of_states = [];

                for ors=1:1:size(orders_of_states, 1)
                   state_order = orders_of_states(ors,:);
                   ordered_states = current_state_assignment(state_order);
                   zero_ordered_states = cell2mat(cellfun(@(x) zero_empty(x), ordered_states, 'UniformOutput', false));
                   perm_id = sum(zero_ordered_states.*((n+1).^(linspace(0, length(zero_ordered_states)-1, length(zero_ordered_states)))));
                   if (any(perms_of_states == perm_id))
                       continue;
                   else
                       perms_of_states = [perms_of_states; perm_id];
                   end
                   if (isempty(ordered_states{r+1}))
                       continue;
                   end
                   % Define state assignment
%                    s = zeros(n, 2*n);
                   s = zeros(n, n);
                   for ors_=1:1:(r+1)
                       s(ordered_pseudo_input{ors_}, ordered_states{ors_}) = 1;
%                        s(ordered_pseudo_input{ors_}, ordered_states{ors_} + n) = 1;
                   end
%                    u_x = [u_x; reshape(p_casc, 1, 2*n), reshape(s, 1, 2*n^2)];
                   u_x = [u_x; reshape(p_casc, 1, 2*n), reshape(s, 1, n^2)];
                end
            end
        end

        % Define action tree for purely decoupled
        p_dec = zeros(n,2);
        for ori_=1:1:(r+1)
            p_dec(current_pseudo_inputs{ori_}, 2) = ori_;
        end

        for sub_s=1:1:size(state_assignments{r}, 1)
            current_state_assignment = state_assignments{r}{sub_s};
            orders_of_states = perms(1:1:(r+1));
            perms_of_states = [];

            for ors=1:1:size(orders_of_states, 1)
               state_order = orders_of_states(ors,:);
               ordered_states = current_state_assignment(state_order);
               zero_ordered_states = cell2mat(cellfun(@(x) zero_empty(x), ordered_states, 'UniformOutput', false));
               perm_id = sum(zero_ordered_states.*((n+1).^(linspace(0, length(zero_ordered_states)-1, length(zero_ordered_states)))));
               if (any(perms_of_states == perm_id))
                   continue;
               else
                   perms_of_states = [perms_of_states; perm_id];
               end
               if any(cellfun(@(x) isempty(x), ordered_states))
                   continue;
               end
               % Define state assignment
%                s = zeros(n, 2*n);
               s = zeros(n, n);
               for ors_=1:1:(r+1)
                   s(ordered_pseudo_input{ors_}, ordered_states{ors_}) = 1;
%                    s(ordered_pseudo_input{ors_}, ordered_states{ors_} + n) = 1;
               end
%                u_x = [u_x; reshape(p_dec, 1, 2*n), reshape(s, 1, 2*n^2)];
               u_x = [u_x; reshape(p_dec, 1, 2*n), reshape(s, 1, n^2)];
            end
        end
    end
end

sys.A = rand(8,8);
sys.B = rand(8,4);
sys.X_DIMS = 8;
sys.U_DIMS = 4;

py_data = load('data/encodings.mat');
py_data_recursive = load('data/encodings_recursive.mat');
u_x = py_data.u_x;

[~, u_u_xa, u_u_xi] = unique(u_x, 'rows', 'stable');
disp(strcat('Number of unique decompositions :', num2str(size(u_u_xa,1))));

action_trees = cell(size(u_x, 1), 1);
ABs = zeros(size(u_x, 1), sys.X_DIMS * (sys.X_DIMS + sys.U_DIMS));
encodings = zeros(size(u_x, 1), sys.U_DIMS * (sys.X_DIMS + sys.U_DIMS));
encoding_adjs = zeros(size(u_x, 1), sys.U_DIMS * (sys.X_DIMS + sys.U_DIMS));
for dd=1:1:size(u_x, 1)
    action_trees{dd} = compute_action_tree(reshape(u_x(dd, 1:2*sys.U_DIMS), sys.U_DIMS,2), ...
                                           reshape(u_x(dd, (1 + 2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS));
    [A_, B_] = pre_process_dyn(sys, action_trees{dd});
    ABs(dd, :) = [reshape(A_, 1, sys.X_DIMS * sys.X_DIMS), reshape(B_, 1, sys.X_DIMS * sys.U_DIMS)];
    encodings(dd, :) = encode_binary(sys, action_trees{dd});
    encoding_adjs(dd, :) = encode_binary_adj(action_trees{dd});
end

[~, u_ABsa, u_ABsi] = unique(ABs, 'rows', 'stable');
[~, u_encodingsa, u_encodingsi] = unique(encodings, 'rows', 'stable');
[~, u_encoding_adjsa, u_encoding_adjsi] = unique(encoding_adjs, 'rows', 'stable');
disp(strcat('Number of unique dynamics :', num2str(size(u_ABsa, 1))));
disp(strcat('Number of unique encodings :', num2str(size(u_encodingsa, 1))));
disp(strcat('Number of unique encoding_adjs :', num2str(size(u_encoding_adjsa, 1))));

mean(abs(py_data.encodings - encoding_adjs), 'all')
mean(abs(py_data_recursive.encodings - encodings), 'all')

model = importKerasNetwork('data/manipulator4dof_elqr_predictor_1000.h5');
encoding_adjs = encoding_adjs';
encoding_adjs = reshape(encoding_adjs, 1, sys.U_DIMS *(sys.X_DIMS + sys.U_DIMS), 1, 4000);
tic;
model_predictions = model.predict(encoding_adjs(:,:,:,1:1000));
time_predict = toc;

%% Functions

function encoding = encode_binary(sys, action_tree)

    X_DIMS = sys.X_DIMS;
    U_DIMS = sys.U_DIMS;
    
    state_dependence = zeros(U_DIMS, X_DIMS);
    dependence_matrix = zeros(U_DIMS, U_DIMS);
    
    for aa =1:1:size(action_tree,1)
        node = action_tree{aa, 1};
        if all(node~=0)
            node_bool = zeros(U_DIMS, 1);
            node_bool(node) = 1;
            state_dependence(logical(node_bool * action_tree{aa, 2})) = 1;
            dependence_matrix(logical(node_bool * node_bool')) = 1;
            
            children = action_tree{aa, end};
            while (~isempty(children))
                child = children{1};
                children = children(2:end,:);
                
                child_id = cellfun(@(x) setequal(x,child), action_tree(:,1));
                child = action_tree(child_id,:);
                child_bool = zeros(U_DIMS, 1);
                child_bool(child{1}) = 1;
                
                state_dependence(logical(node_bool * child{2})) = 1;
                dependence_matrix(logical(child_bool * node_bool')) = 1;
                
                if (~isempty(child{end}))
                    children = cat(1, children, child{end});
                end
            end
        end
    end
    
    dependence_matrix = dependence_matrix - diag(diag(dependence_matrix));
    encoding = [reshape(state_dependence, 1, U_DIMS * X_DIMS), ...
                reshape(dependence_matrix, 1, U_DIMS * U_DIMS)];
end

function action_tree = compute_action_tree(p, s)
    
    assert(size(p,1) == size(s,1), 'Not a valid decomposition');

    % Build action tree
    X_DIMS = size(s, 2);
    U_DIMS = size(s, 1);
    
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
end

function [A_, B_] = pre_process_dyn(sys, action_tree)
    
    A = sys.A;
    B = sys.B;
    X_DIMS = sys.X_DIMS;
    U_DIMS = sys.U_DIMS;
    
    A_ = zeros(X_DIMS, X_DIMS);
    B_ = zeros(X_DIMS, U_DIMS);
    
    sub_states = cellfun(@(x) zeros(1, X_DIMS), action_tree(:,1), 'UniformOutput', false);
    sub_actions = cellfun(@(x) zeros(1, U_DIMS), action_tree(:,1), 'UniformOutput', false);
    action_tree = cat(2, action_tree(:, 1:2), sub_states, sub_actions, action_tree(:,3:end));
    
    while (size(action_tree, 1) > 1)
        % Find leaf nodes
        leaf_node_ids = cellfun(@(x) isempty(x), action_tree(:,end));
        leaf_nodes = action_tree(leaf_node_ids, :);
        
        for ii=1:1:size(leaf_nodes,1)
            leaf_node = leaf_nodes{ii,1};
            leaf_states = leaf_nodes{ii,2};
            leaf_substates = leaf_nodes{ii,3};
            leaf_subactions = leaf_nodes{ii,4};
            
%             if (isfield(sys, 'fxfu_func'))
%                 fxfu = sys.fxfu_func(sys.l_point, u0_);
%             else
%                 fxfu = eval(subs(sys.fxfu, [sys.xu], [sys.l_point; u0_]));
%             end
            
            state_mask = logical(leaf_states' * leaf_states);
            A_(state_mask) = A(state_mask);
            substate_mask = logical(leaf_states' * leaf_substates);
            A_(substate_mask) = A(substate_mask);
            
            leaf_node_bool = zeros(1, U_DIMS);
            leaf_node_bool(leaf_node) = 1;
            state_action_mask = logical(leaf_states' * leaf_node_bool);
            B_(state_action_mask) = B(state_action_mask);
            state_subaction_mask = logical(leaf_states' * leaf_subactions);
            B_(state_subaction_mask) = B(state_subaction_mask);
            
            % Find parent and update action tree
            parent_input = leaf_nodes{ii, 5};
            parent_node_id = find(cellfun(@(x) isempty(setdiff(x, parent_input)) && isempty(setdiff(parent_input, x)),...
                                  action_tree(:, 1)));
            assert(length(parent_node_id)==1, 'Invalid Parent Node ID');
            parent_node = action_tree(parent_node_id, :);

            assert(all(~(parent_node{2} .* (leaf_states + leaf_substates))), 'Invalid state overlap');
            assert(all(~(parent_node{3} .* (leaf_states + leaf_substates))), 'Invalid state substate overlap');
            parent_node{3} = parent_node{3} + leaf_states + leaf_substates;
            parent_node{4} = parent_node{4} + leaf_node_bool + leaf_subactions;

            % Find and delete the leaf_node in the list of children
            children_list = parent_node{end};
            childID = cellfun(@(x) isempty(setdiff(x, leaf_node)) ...
                                   && isempty(setdiff(leaf_node, x)), ...
                            children_list);
            children_list(childID) = [];
            parent_node{end} = children_list;

            action_tree(parent_node_id, :) = parent_node;
        end
        action_tree(leaf_node_ids, :) = [];
    end
end

function x = zero_empty(x)
    if (isempty(x))
       x = 0;
    end
end