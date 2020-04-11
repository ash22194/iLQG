function ddl = ddldudx_CartPoleFirst(sys, x, u)
    
    X_DIMS_FREE = sys.X_DIMS_FREE;
    U_DIMS_FREE = sys.U_DIMS_FREE;
    ddl = zeros(length(U_DIMS_FREE), ...
                length(X_DIMS_FREE), ...
                size(x, 2)); 
end