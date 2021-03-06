function cost = l_Biped2DSecond(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
%     U_DIMS_FIXED = linspace(1,2,2)';
%     U_DIMS_FIXED(U_DIMS_FREE) = [];
    lims = sys.lims;
    action = zeros(length(U_DIMS_FREE) + length(U_DIMS_FIXED), size(u, 2));
    action(U_DIMS_FREE, :) = u;
    for ii=1:1:size(u, 2)
        action(U_DIMS_FIXED, ii) = k(:, ii) + K(:,:, ii)*(x(:, ii) - xn(:, ii));
        action(:,ii) = max(lims(:,1), min(lims(:,2), action(:, ii)));
    end

    goal = [sys.goal(1); sys.goal(2); sys.goal(3); sys.goal(4); sys.goal(5); sys.goal(6)];
    cost = diag((x - goal)'*sys.Q*(x - goal) + action'*sys.R*action)';
    
end