function dl = dldx_CartPoleFirst_WrapAround(sys, x, u)
    
    X_DIMS_FREE = sys.X_DIMS_FREE;
    wraparound_dim = find(sys.X_DIMS_FREE==3);
    x(wraparound_dim,:) = mod(x(wraparound_dim,:), 2*pi);
    
    goal = sys.goal(X_DIMS_FREE);
    Q = sys.Q(X_DIMS_FREE, :);
    Q = Q(:, X_DIMS_FREE);

    dl = 2 * Q * (x - goal);
    
end