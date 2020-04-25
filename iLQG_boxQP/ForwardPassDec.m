function [X, U] = ForwardPassDec(k, K, Xn, U_DIMS, x_start, dynamics_discrete)
NUM_CTRL = size(k, 2)-1;
U_DIM = size(k, 1);
X_DIM = size(x_start, 1);

X = zeros(X_DIM, NUM_CTRL+1);
U = ones(U_DIM, NUM_CTRL);

X(:, 1) = x_start;
for t = 1 : NUM_CTRL
    U(:, t) = zeros(U_DIM, 1); 
    for ii=1:1:size(U_DIMS, 1)
        U_DIMS_ii = U_DIMS{ii};
        [~, t_closest] = min(vecnorm(X(:, t) - Xn(((ii-1)*X_DIM+1):ii*X_DIM, :), 2, 1));
        U(U_DIMS_ii, t) = k(U_DIMS_ii, t_closest)...
                   + K(U_DIMS_ii, :, t_closest)*(X(:, t) - Xn(((ii-1)*X_DIM+1):ii*X_DIM, t_closest));
    end
    [X(:, t+1), ~] = dynamics_discrete(X(:, t), U(:, t));
end
end
