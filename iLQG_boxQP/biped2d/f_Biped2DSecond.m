function dx = f_Biped2DSecond(sys, x, u, k, K, xn)

    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
    
    lims = sys.lims;
    U = zeros(4, size(u,2));
    U(U_DIMS_FREE, :) = u;

    TRAJ_LENGTH = size(x,2);
    x_ = x - xn;
    x_ = reshape(x_, 1, size(x_,1), TRAJ_LENGTH);
    x_ = repmat(x_, [length(U_DIMS_FIXED), 1, 1]);
    K_ = reshape(sum(K.*x_, 2), length(U_DIMS_FIXED), TRAJ_LENGTH);
    U(U_DIMS_FIXED, :) = k + K_(:,1,:);
    U = max(repmat(lims(:,1), [1, size(u,2)]), ...
            min(repmat(lims(:,2), [1, size(u,2)]), U));

    ca1       = cos(x(2,:));
    sa1       = sin(x(2,:));
    x_hip     = x(1,:).*ca1;
    z_hip     = x(1,:).*sa1;
    l2        = sqrt((x_hip + sys.df).^2 + z_hip.^2);
    a2        = acos((sys.df + x_hip)./l2);
    ca2       = cos(a2);
    sa2       = sin(a2);
    contact1  = (x(1,:)<=sys.l0);
    contact2  = (l2<=sys.l0);
    F1  = U(1,:).*contact1;
    F2  = U(2,:).*contact2;
    Fo1 = U(3,:).*contact1./x(1,:);
    Fo2 = U(4,:).*contact2./l2;
    
    x_hip_dot = x(3,:) + sys.d*cos(x(5,:)).*x(6,:);
    z_hip_dot = x(4,:) + sys.d*sin(x(5,:)).*x(6,:);
    dx = [x_hip_dot.*ca1 + z_hip_dot.*sa1;
          (-x_hip_dot.*sa1 + z_hip_dot.*ca1)./x(1,:);
          (F1.*ca1 + Fo1.*sa1 + F2.*ca2 + Fo2.*sa2)/sys.m;
          ((F1.*sa1 - Fo1.*ca1 + F2.*sa2 - Fo2.*ca2)/sys.m - sys.g);
          x(6,:);
          (Fo1.*(x(1,:) + sys.d*sin(x(2,:) - x(5,:))) ...
             + Fo2.*(l2 + sys.d*sin(a2 - x(5,:))) ...
             + F1.*(sys.d*cos(x(2,:) - x(5,:))) ...
             + F2.*(sys.d*cos(a2 - x(5,:))))/sys.I];

end
