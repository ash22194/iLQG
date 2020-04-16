function cost = l_CartPoleFirst_WrapAround(sys, x, u)
    
    X_DIMS_FREE = sys.X_DIMS_FREE;
    U_DIMS_FREE = sys.U_DIMS_FREE;
    wraparound_dim = find(sys.X_DIMS_FREE==3);
    x(wraparound_dim,:) = mod(x(wraparound_dim,:), 2*pi);
    
    goal = sys.goal(X_DIMS_FREE);
%     lims = sys.lims(U_DIMS_FREE, :);
%     action = max(lims(:,1)*ones(1,size(u,2)), ...
%                  min(lims(:,2)*ones(1,size(u,2)), u));
    action = u;
    Q = sys.Q(X_DIMS_FREE, :);
    Q = Q(:, X_DIMS_FREE);
    R = sys.R(U_DIMS_FREE, :);
    R = R(:, U_DIMS_FREE);
    cost = diag((x - goal)'*Q*(x - goal) + action'*R*action)';
    
end