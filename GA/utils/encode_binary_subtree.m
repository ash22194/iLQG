function encoding = encode_binary_subtree(p, s)

U_DIMS = size(s, 1);
X_DIMS = size(s, 2);

adj_matrix = false(U_DIMS);
for aa =1:1:U_DIMS
    node_bool = all(p(aa, :) == p, 2);
    adj_matrix(node_bool * node_bool') = true;

    child_bool = any(p(:,1) == find(node_bool)', 2);
    adj_matrix(child_bool * node_bool') = true;
end
adj_matrix = adj_matrix - diag(diag(adj_matrix));
adjacency_matrix = false(U_DIMS);

s = logical(s);

nodes = 1:1:U_DIMS;
leafnodebool = ~any(adj_matrix, 1);
num_leaf_nodes = sum(leafnodebool);
while(num_leaf_nodes > 0)
    leafnode = nodes(leafnodebool);
    for ll=1:1:length(leafnode)
        parent = adj_matrix(leafnode(ll), :);
        if (sum(parent) > 0)
            parent = nodes(parent);
            adjacency_matrix(parent, parent) = true;
            adjacency_matrix(parent, :) ...
               = adjacency_matrix(parent, :) ...
                 | repmat(adjacency_matrix(leafnode(ll), :), [length(parent), 1]);
            s(parent, :) = s(parent, :) | repmat(s(leafnode(ll), :), [length(parent), 1]);
        end
    end
    
    nodes(leafnodebool) = [];
    adj_matrix = adj_matrix(:, ~leafnodebool);
    leafnodebool = ~any(adj_matrix, 1);
    num_leaf_nodes = sum(leafnodebool);
end

encoding = [reshape(s, 1, U_DIMS*X_DIMS), reshape(adjacency_matrix, 1, U_DIMS^2)];
end