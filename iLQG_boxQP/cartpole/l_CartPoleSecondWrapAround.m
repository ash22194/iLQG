function cost = l_CartPoleSecondWrapAround(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = linspace(1,2,2)';
    U_DIMS_FIXED(U_DIMS_FREE) = [];
    lims = sys.lims;
    action = zeros(2, size(u, 2));
    action(U_DIMS_FREE, :) = u;
    for ii=1:1:size(u, 2)
        action(U_DIMS_FIXED, ii) = k(:, ii) + K(:,:, ii)*(x(:, ii) - xn(:, ii));
        action(:,ii) = max(lims(:,1), min(lims(:,2), action(:, ii)));
    end
    
    while(any(x(3, :) > 2*pi))
        greaterThan = x(3,:) > 2*pi;
        x(3,greaterThan) = x(3,greaterThan) - 2*pi;
    end
    while(any(x(3, :) < 0))
        lessThan = x(3,:) < 0;
        x(3,lessThan) = x(3,lessThan) + 2*pi;
    end
    
    goal = [sys.goal(1); sys.goal(2); sys.goal(3); sys.goal(4)];
    cost = diag((x - goal)'*sys.Q*(x - goal) + action'*sys.R*action)';
    
end