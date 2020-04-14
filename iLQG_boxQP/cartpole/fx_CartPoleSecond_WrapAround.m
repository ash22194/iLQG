function fx = fx_CartPoleSecond_WrapAround(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
%     U_DIMS_FIXED = linspace(1, 2, 2)';
%     U_DIMS_FIXED(U_DIMS_FREE) = [];
    fx = zeros(4, 4, size(x, 2));
    lims = sys.lims;
    while (any(x(3,:) > 2*pi))
        gthan = x(3,:) > 2*pi;
        x(3,gthan) = x(3,gthan) - 2*pi;
    end
    while (any(x(3,:) < 0))
        lthan = x(3,:) < 0;
        x(3,lthan) = x(3,lthan) + 2*pi;
    end
    
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