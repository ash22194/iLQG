function dl = dldx_CartPoleFirst_WrapAround(sys, x, u)
    
    X_DIMS_FREE = sys.X_DIMS_FREE;
    wraparound_dim = find(sys.X_DIMS_FREE==3);
    while (any(x(wraparound_dim,:) > 2*pi))
        gthan = x(wraparound_dim,:) > 2*pi;
        x(wraparound_dim,gthan) = x(wraparound_dim,gthan) - 2*pi;
    end
    while (any(x(wraparound_dim,:) < 0))
        lthan = x(wraparound_dim,:) < 0;
        x(wraparound_dim,lthan) = x(wraparound_dim,lthan) + 2*pi;
    end
    
    goal = sys.goal(X_DIMS_FREE);
    Q = sys.Q(X_DIMS_FREE, :);
    Q = Q(:, X_DIMS_FREE);

    dl = 2 * Q * (x - goal);
    
end