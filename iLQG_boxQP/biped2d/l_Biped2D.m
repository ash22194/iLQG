function cost = l_Biped2D(sys, x, u)
    
    pos = [x(1,:); x(2,:); x(3,:); x(4,:); x(5,:); x(6,:)];
    goal = [sys.goal(1); sys.goal(2); sys.goal(3); sys.goal(4); sys.goal(5); sys.goal(6)];
    action = [u(1,:); u(2,:); u(3,:); u(4,:)];
    cost = diag((pos - goal)'*sys.Q*(pos - goal) + action'*sys.R*action)';
    
end