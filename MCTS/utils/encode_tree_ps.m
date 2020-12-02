function encoding = encode_tree_ps(action_tree)
    
    X_DIMS = sum(cellfun(@(x) sum(x), action_tree(:,2)));
    U_DIMS = sum(cellfun(@(x) length(x), action_tree(:,1))) - 1;
    state_assignment = zeros(U_DIMS, X_DIMS);
    input_coupling = zeros(U_DIMS, U_DIMS);
    input_dependence = zeros(U_DIMS, U_DIMS);
    
    % action tree = {input_tuple, state_dims, }
    for ii=1:1:size(action_tree, 1)
        if (all(action_tree{ii, 1}~=0))
            state_assignment(action_tree{ii,1}, :) = ones(length(action_tree{ii,1}), 1) ...
                                                      * double(action_tree{ii, 2});

            input_coupling(action_tree{ii, 1}, action_tree{ii, 1}) = 1;

            children = set_union(action_tree{ii, end}{:});
            input_dependence(action_tree{ii, 1}, children) = 1;
        end
    end
    
    encoding = [reshape(input_coupling, 1, U_DIMS^2), ...
                reshape(input_dependence, 1, U_DIMS^2), ...
                reshape(state_assignment, 1, U_DIMS*X_DIMS)];
end

function final_set = set_union(varargin)
    
    final_set = [];
    for ii=1:nargin
        final_set = union(final_set, varargin{ii});
    end
end