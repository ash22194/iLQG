function ddl = ddldx_Biped2DFirst(sys, x, u)
    
    X_DIMS_FREE = sys.X_DIMS_FREE;
    Q = sys.Q(X_DIMS_FREE, :);
    Q = Q(:, X_DIMS_FREE);
    ddl = repmat(2 * Q, [1,1,size(x, 2)]);

end