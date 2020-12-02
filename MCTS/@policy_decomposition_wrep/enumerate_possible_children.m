function all_possible_children = enumerate_possible_children(de)
            
    X_DIMS = de.sys.X_DIMS;
    U_DIMS = de.sys.U_DIMS;
    all_possible_children = cell(de.childnode_maxcount, 1);
    children_count = 0;
    for jj=1:1:size(de.action_tree,1)
        inputs = de.action_tree{jj, 1};
        num_inputs = length(inputs);
        if (num_inputs==1)
            continue;
        end
        % Possible non-empty 2-groupings of the inputs
        input_2groups = cell(round(distr_n_into_r(num_inputs, 2) / 2), 2);
        num_input_groups = 0;
        for ii=1:1:floor(num_inputs/2)
            size_ii_groups = nchoosek(inputs, ii);
            if (ii == ceil(num_inputs/2))
                size_ii_groups = size_ii_groups(1:(size(size_ii_groups, 1)/2),:);
            end
            num_ii_groups  = size(size_ii_groups, 1);
            size_ii_groups = cellfun(@(x) x', mat2cell(size_ii_groups, ones(num_ii_groups, 1), ii), 'UniformOutput', false);

            input_2groups((num_input_groups+1):(num_input_groups+num_ii_groups), 1) = size_ii_groups;
            input_2groups((num_input_groups+1):(num_input_groups+num_ii_groups), 2) ...
                = cellfun(@(x) setdiff(inputs, x), size_ii_groups, 'UniformOutput', false);

            num_input_groups = num_input_groups + num_ii_groups;
        end

        states = de.action_tree{jj, 2};
        num_states = sum(states);
        states = find(states);
        % Possible 2-groupings of the states
        state_2groups = cell((2^(num_states)), 2);
        num_state_groups = 0;
        for ss=0:1:num_states
            if (num_states ~= 1)
                size_ss_groups = nchoosek(states, ss);
            else
                if (ss==0)
                    size_ss_groups = zeros(1,0);
                else
                    size_ss_groups = states;
                end
            end
            num_ss_groups  = size(size_ss_groups, 1);
            size_ss_groups = mat2cell(size_ss_groups, ones(num_ss_groups, 1), ss);

            state_2groups((num_state_groups+1):(num_state_groups+num_ss_groups), 1) = size_ss_groups;
            state_2groups((num_state_groups+1):(num_state_groups+num_ss_groups), 2) ...
                = cellfun(@(x) setdiff(states, x), size_ss_groups, 'UniformOutput', false);

            num_state_groups = num_state_groups + num_ss_groups;
        end
        
        num_children = length(de.action_tree{jj, end});
        children = linspace(1, num_children, num_children);
        % Possible 2-groupings of the children
        children_2groups = cell((2^(num_children)), 2);
        num_children_groups = 0;
        for ch=0:1:num_children
            if (num_children ~= 1)
                size_ch_groups = nchoosek(children, ch);
            else
                if (ch==0)
                    size_ch_groups = zeros(1,0);
                else
                    size_ch_groups = ones(1,1);
                end
            end
            num_ch_groups  = size(size_ch_groups, 1);
            size_ch_groups = cellfun(@(x) x', mat2cell(size_ch_groups, ones(num_ch_groups, 1), ch), 'UniformOutput', false);
            
            children_2groups((num_children_groups+1):(num_children_groups+num_ch_groups), 1) ...
                = cellfun(@(x) de.action_tree{jj, end}(x), size_ch_groups, 'UniformOutput', false);
            children_2groups((num_children_groups+1):(num_children_groups+num_ch_groups), 2) ...
                = cellfun(@(x) de.action_tree{jj, end}(setdiff(children, x)'), size_ch_groups, 'UniformOutput', false);

            num_children_groups = num_children_groups + num_ch_groups;
        end
        
        p__ = de.p;
        s__ = de.s;
        parent_node = p__(inputs, 1);
        assert(all(parent_node(1)==parent_node, 'all'), 'Coupled inputs must have the same parent!');
        parent_node = parent_node(1);
        
        curr_child_id = p__(inputs, 2);
        assert(all(curr_child_id(1)==curr_child_id, 'all'), 'Coupled inputs must have the same child ID!');
        curr_child_id = curr_child_id(1);
        all_child_ids = p__(p__(:,1) == parent_node, 2);
        
        for ni=1:1:num_input_groups
            curr_input_group = input_2groups(ni, :); % Cell array of size 2
            
            for nc=1:1:num_children_groups
                % Distribute the child nodes
                curr_children_group = children_2groups(nc, :); % Cell array of size 2
                num_ch1 = size(curr_children_group{1},1);
                num_ch2 = size(curr_children_group{2},1);
                
                for ns=1:1:num_state_groups
                    curr_state_group = state_2groups(ns, :);

                    if ~(((num_ch1==0) && (isempty(curr_state_group{1}))) ...
                        || ((num_ch2==0) && (isempty(curr_state_group{2}))))
                        % Decoupled
                        % Both input groups have same parent but different child IDs
                        p_ = p__;
                        p_(curr_input_group{1}, 2) = curr_child_id;
                        p_(curr_input_group{2}, 2) = max(all_child_ids) + 1;
                        for ncc=1:1:length(curr_children_group{1})
                            p_(curr_children_group{1}{ncc}, 1) = curr_input_group{1}(1);
                            p_(curr_children_group{1}{ncc}, 2) = ncc;
                        end
                        for ncc=1:1:length(curr_children_group{2})
                            p_(curr_children_group{2}{ncc}, 1) = curr_input_group{2}(1);
                            p_(curr_children_group{2}{ncc}, 2) = ncc;
                        end
                        
                        s_ = s__;
                        s_(curr_input_group{1}, :) = false;
                        s_(curr_input_group{1}, curr_state_group{1}) = true;
                        s_(curr_input_group{2}, :) = false;
                        s_(curr_input_group{2}, curr_state_group{2}) = true;
                        
                        all_possible_children{children_count+1} = [reshape(p_, 1, 2*U_DIMS),...
                                                                   reshape(s_, 1, U_DIMS*X_DIMS)];
                        children_count = children_count + 1;
                    end

                    if ~((num_ch2==0) && isempty(curr_state_group{2}))
                        % Cascaded 1 after 2
                        % Input group 1 is parent for input group 2
                        p_ = p__;
                        p_(curr_input_group{2}, 1) = curr_input_group{1}(1);
                        p_(curr_input_group{1}, 2) = curr_child_id;
                        for ncc=1:1:length(curr_children_group{1})
                            p_(curr_children_group{1}{ncc}, 1) = curr_input_group{1}(1);
                            p_(curr_children_group{1}{ncc}, 2) = ncc;
                        end
                        p_(curr_input_group{2}, 2) = length(curr_children_group{1}) + 1;
                        for ncc=1:1:length(curr_children_group{2})
                            p_(curr_children_group{2}{ncc}, 1) = curr_input_group{2}(1);
                            p_(curr_children_group{2}{ncc}, 2) = ncc;
                        end
                        
                        s_ = s__;
                        s_(curr_input_group{1}, :) = false;
                        s_(curr_input_group{1}, curr_state_group{1}) = true;
                        s_(curr_input_group{2}, :) = false;
                        s_(curr_input_group{2}, curr_state_group{2}) = true;
                        
                        all_possible_children{children_count+1} = [reshape(p_, 1, 2*U_DIMS),...
                                                                   reshape(s_, 1, U_DIMS*X_DIMS)];
                        children_count = children_count + 1;
                    end

                    if ~((num_ch1==0) && isempty(curr_state_group{1}))
                        % Cascaded 2 after 1
                        % Input group 2 is parent for input group 1
                        p_ = p__;
                        p_(curr_input_group{1}, 1) = curr_input_group{2}(1);
                        p_(curr_input_group{2}, 2) = curr_child_id;
                        for ncc=1:1:length(curr_children_group{2})
                            p_(curr_children_group{2}{ncc}, 1) = curr_input_group{2}(1);
                            p_(curr_children_group{2}{ncc}, 2) = ncc;
                        end
                        p_(curr_input_group{1}, 2) = length(curr_children_group{2}) + 1;
                        for ncc=1:1:length(curr_children_group{1})
                            p_(curr_children_group{1}{ncc}, 1) = curr_input_group{1}(1);
                            p_(curr_children_group{1}{ncc}, 2) = ncc;
                        end
                        
                        s_ = s__;
                        s_(curr_input_group{1}, :) = false;
                        s_(curr_input_group{1}, curr_state_group{1}) = true;
                        s_(curr_input_group{2}, :) = false;
                        s_(curr_input_group{2}, curr_state_group{2}) = true;
                        
                        all_possible_children{children_count+1} = [reshape(p_, 1, 2*U_DIMS),...
                                                                   reshape(s_, 1, U_DIMS*X_DIMS)];
                        children_count = children_count + 1;
                    end
                end
            end
        end
    end
    
    assert(children_count==de.childnode_maxcount, 'Mismatch in the number of expected and enumerated child nodes');
    de.childnodes = all_possible_children;
end