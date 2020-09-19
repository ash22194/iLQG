function c = cost(sys, x, u)
    %% Cost calculations
    
    c = diag((x - sys.l_point)'*sys.Q*(x - sys.l_point) ...
             + (u - sys.u0)'*sys.R*(u - sys.u0))';

end