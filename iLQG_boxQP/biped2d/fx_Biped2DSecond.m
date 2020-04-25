function fx = fx_Biped2DSecond(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
    U_DIM = length(U_DIMS_FREE) + length(U_DIMS_FIXED);
    
    X_DIMS_FREE = sys.X_DIMS_FREE;
    X_DIMS_FIXED = sys.X_DIMS_FIXED;
    X_DIM = length(X_DIMS_FREE) + length(X_DIMS_FIXED);

    fx = zeros(X_DIM, X_DIM, size(x, 2)); % TODO: What if there are more than 2 layers of cascade?
    lims = sys.lims;
    for ii=1:1:size(u, 2)
        U = zeros(U_DIM, 1);
        U(U_DIMS_FREE, 1) = u(:, ii);
        U(U_DIMS_FIXED, 1) = k(:, ii) + K(:,:, ii) * (x(:, ii) - xn(:, ii));
        U = max(lims(:,1), min(lims(:,2), U));
        
        K_ = zeros(U_DIM, X_DIM);
        K_(U_DIMS_FIXED, :) = K(:,:, ii);
        fx(:,:, ii) = fx_Biped2D(sys, x(:, ii), U) ...
                      + fu_Biped2D(sys, x(:, ii), U) * K_;
    end
end