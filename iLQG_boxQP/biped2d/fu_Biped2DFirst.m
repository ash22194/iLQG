function fu = fu_Biped2DFirst(sys, x, u)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
    l_point = sys.l_point;
    X_DIMS_FREE = sys.X_DIMS_FREE;
    X_DIMS_FIXED = sys.X_DIMS_FIXED;

    X = zeros(length(X_DIMS_FREE) +  length(X_DIMS_FIXED), size(x, 2));
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
    ca1       = ca1.*contact1;
    sa1       = sa1.*contact1;
    contact2  = (l2<=sys.l0);
    ca2       = ca2.*contact2;
    sa2       = sa2.*contact2;
    
    fu = zeros(size(X, 1), length(U_DIMS_FREE) +  length(U_DIMS_FIXED), size(X, 2));
    fu(3,1,:) = ca1/sys.m;
    fu(3,2,:) = ca2/sys.m;
    fu(3,3,:) = sa1./(X(1,:)*sys.m);
    fu(3,4,:) = sa2./(l2*sys.m);
    
    fu(4,1,:) = sa1/sys.m;
    fu(4,2,:) = sa2/sys.m;
    fu(4,3,:) = -ca1./(X(1,:)*sys.m);
    fu(4,4,:) = -ca2./(l2*sys.m);
    
    fu(6,1,:) = contact1.*(sys.d*cos(X(2,:) - X(5,:)))/sys.I;
    fu(6,2,:) = contact2.*(sys.d*cos(a2 - X(5,:)))/sys.I;
    fu(6,3,:) = contact1.*(1 + sys.d./X(1,:).*sin(X(2,:) - X(5,:)))/sys.I;
    fu(6,4,:) = contact2.*(1 + sys.d./l2.*sin(a2 - X(5,:)))/sys.I;

    fu = fu(X_DIMS_FREE, :, :);
    fu = fu(:, U_DIMS_FREE, :);
end
