function dx = cartpole_dyn_gridbased(t, x, policyF, policyT, limits, sys)
    
    x_ = [min(limits(1, 2), max(limits(1, 1), x(1)));
          min(limits(2, 2), max(limits(2, 1), x(2)));
          min(limits(3, 2), max(limits(3, 1), x(3)));
          min(limits(4, 2), max(limits(4, 1), x(4)))];
    

    u_pole = policyT(x_(1), x_(2), x_(3), x_(4));
    u_cart = policyF(x_(1), x_(2), x_(3), x_(4));

    dx = zeros(4,1);
    dx(1) = x(2);
    dx(2) = (u_cart - u_pole*cos(x(3))/sys.l...
            + sys.mp*sys.l*(x(4)^2)*sin(x(3)) + 0.5*sys.mp*sys.g*sin(2*x(3)))/(sys.mc + sys.mp*sin(x(3))^2);
    dx(3) = x(4);
    dx(4) = (u_pole/sys.l*(sys.mc/sys.mp+1) - u_cart*cos(x(3))...
             - 0.5*sys.mp*sys.l*(x(4)^2)*sin(2*x(3)) - (sys.mc+sys.mp)*sys.g*sin(x(3)))/(sys.l*(sys.mc + sys.mp*sin(x(3))^2));

end