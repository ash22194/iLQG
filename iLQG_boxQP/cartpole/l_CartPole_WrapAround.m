function cost = l_CartPole_WrapAround(sys, x, u)
    x(3,:) = mod(x(3,:), 2*pi);
    pos = [x(1, :); x(2, :); x(3, :); x(4, :)];
    goal = [sys.goal(1); sys.goal(2); sys.goal(3); sys.goal(4)];
    action = [u(1, :); u(2, :)];
    cost = diag((pos - goal)'*sys.Q*(pos - goal) + action'*sys.R*action)';
    
end