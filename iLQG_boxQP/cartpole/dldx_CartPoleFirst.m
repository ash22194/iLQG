function dl = dldx_CartPoleFirst(sys, x, u)
    
    X_DIMS_FREE = sys.X_DIMS_FREE;
    
    goal = sys.goal(X_DIMS_FREE);
    Q = sys.Q(X_DIMS_FREE, :);
    Q = Q(:, X_DIMS_FREE);

    dl = 2 * Q * (x - goal);
    
end