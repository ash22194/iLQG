function dx = f_Biped2DFirst(sys, x, u)

    l_point = sys.l_point;

    X_DIMS_FREE = sys.X_DIMS_FREE;
    X_DIMS_FIXED = sys.X_DIMS_FIXED;

    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;

    U = zeros(length(U_DIMS_FIXED) + length(U_DIMS_FREE), size(u, 2));
    U(U_DIMS_FREE, :) = u;
    X = zeros(length(X_DIMS_FIXED) + length(X_DIMS_FREE), size(x, 2));
    X(X_DIMS_FREE, :) = x;
    X(X_DIMS_FIXED, :) = repmat(l_point(X_DIMS_FIXED), [1, size(x,2)]);

    ca1       = cos(X(2,:));
    sa1       = sin(X(2,:));
    x_hip     = X(1,:).*ca1;
    z_hip     = X(1,:).*sa1;
    l2        = sqrt((x_hip + sys.df).^2 + z_hip.^2);
    a2        = acos((sys.df + x_hip)./l2);
    ca2       = cos(a2);
    sa2       = sin(a2);
    contact1  = (X(1,:)<=sys.l0);
    contact2  = (l2<=sys.l0);
    F1  = U(1,:).*contact1;
    F2  = U(2,:).*contact2;
    Fo1 = U(3,:).*contact1./x(1,:);
    Fo2 = U(4,:).*contact2./l2;
    
    x_hip_dot = X(3,:) + sys.d*cos(X(5,:)).*X(6,:);
    z_hip_dot = X(4,:) + sys.d*sin(X(5,:)).*X(6,:);
    dx = [x_hip_dot.*ca1 + z_hip_dot.*sa1;
          (-x_hip_dot.*sa1 + z_hip_dot.*ca1)./X(1,:);
          (F1.*ca1 + Fo1.*sa1 + F2.*ca2 + Fo2.*sa2)/sys.m;
          ((F1.*sa1 - Fo1.*ca1 + F2.*sa2 - Fo2.*ca2)/sys.m - sys.g);
          X(6,:);
          (Fo1.*(X(1,:) + sys.d*sin(X(2,:) - X(5,:))) ...
             + Fo2.*(l2 + sys.d*sin(a2 - X(5,:))) ...
             + F1.*(sys.d*cos(X(2,:) - X(5,:))) ...
             + F2.*(sys.d*cos(a2 - X(5,:))))/sys.I];
    
    dx = dx(X_DIMS_FREE);
end
