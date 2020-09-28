function c = cost(sys, x, u)
    %% Cost calculations
    
    x_bar = sign(x - sys.l_point).*mod(abs(x - sys.l_point), sys.cxmod);
    c = diag(x_bar'*sys.Q*x_bar ...
             + (u - sys.u0)'*sys.R*(u - sys.u0))';

end