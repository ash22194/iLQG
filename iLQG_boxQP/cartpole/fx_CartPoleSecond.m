function fx = fx_CartPoleSecond(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
%     U_DIMS_FIXED = linspace(1, 2, 2)';
%     U_DIMS_FIXED(U_DIMS_FREE) = [];
    fx = zeros(4, 4, size(x, 2));
    lims = sys.lims;
    for ii=1:1:size(u, 2)
        U = zeros(2, 1);
        U(U_DIMS_FREE, 1) = u(:, ii);
        U(U_DIMS_FIXED, 1) = k(:, ii) + K(:,:, ii) * (x(:, ii) - xn(:, ii));
        U = max(lims(:,1), min(lims(:,2), U));
        
        K_ = zeros(2, 4);
        K_(U_DIMS_FIXED, :) = K(:,:, ii);
        fx(:,:, ii) = fx_CartPole(sys, x(:, ii), U) ...
                      + fu_CartPole(sys, x(:, ii), U) * K_;
    end
end