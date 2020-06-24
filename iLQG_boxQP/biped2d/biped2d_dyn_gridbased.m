function dx = biped2d_dyn_gridbased(t, x, policy1, policy2, policy3, policy4, limits, sys)
    
    x_ = [min(limits(1, 2), max(limits(1, 1), x(1)));
          min(limits(2, 2), max(limits(2, 1), x(2)));
          min(limits(3, 2), max(limits(3, 1), x(3)));
          min(limits(4, 2), max(limits(4, 1), x(4)));
          min(limits(5, 2), max(limits(5, 1), x(5)));
          min(limits(6, 2), max(limits(6, 1), x(6)))];
    
    F1 = policy1(x_(1), x_(2), x_(3), x_(4), x_(5), x_(6));
    F2 = policy2(x_(1), x_(2), x_(3), x_(4), x_(5), x_(6));
    T1 = policy3(x_(1), x_(2), x_(3), x_(4), x_(5), x_(6));
    T2 = policy4(x_(1), x_(2), x_(3), x_(4), x_(5), x_(6));
    
    ca1       = cos(x(2));
    sa1       = sin(x(2));
    x_hip     = x(1)*ca1;
    z_hip     = x(1)*sa1;
    l2        = sqrt((x_hip + sys.df)^2 + z_hip^2);
    a2        = acos((sys.df + x_hip)/l2);
    ca2       = cos(a2);
    sa2       = sin(a2);
    contact1  = (x(1)<=sys.l0);
    contact2  = (l2<=sys.l0);
    F1 = F1*contact1;
    Fo1 = T1*contact1/x(1);
    F2 = F2*contact2;
    Fo2 = T2*contact2/l2;
    
    dx = zeros(6,1);
    x_hip_dot = x(3) + sys.d*cos(x(5))*x(6);
    z_hip_dot = x(4) + sys.d*sin(x(5))*x(6);
    dx(1) = x_hip_dot*ca1 + z_hip_dot*sa1;
    dx(2) = (-x_hip_dot*sa1 + z_hip_dot*ca1)/x(1);
    dx(3) = (F1*ca1 + Fo1*sa1 ...
               + F2*ca2 + Fo2*sa2)/sys.m;
    dx(4) = ((F1*sa1 - Fo1*ca1 ...
               + F2*sa2 - Fo2*ca2)/sys.m - sys.g);
    dx(5) = x(6);
    dx(6) = (Fo1*(x(1) + sys.d*sin(x(2) - x(5))) ...
             + Fo2*(l2 + sys.d*sin(a2 - x(5))) ...
             + F1*sys.d*cos(x(2) - x(5)) ...
             + F2*sys.d*cos(a2 - x(5)))/sys.I;

end