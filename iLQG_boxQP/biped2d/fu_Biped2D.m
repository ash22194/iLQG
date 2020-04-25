function fu = fu_Biped2D(sys, x, u)
    
    ca1       = cos(x(2,:));
    sa1       = sin(x(2,:));
    x_hip     = x(1,:).*ca1;
    z_hip     = x(1,:).*sa1;
    l2        = sqrt((x_hip + sys.df).^2 + z_hip.^2);
    a2        = acos((sys.df + x_hip)./l2);
    ca2       = cos(a2);
    sa2       = sin(a2);
    contact1  = (x(1,:)<=sys.l0);
    ca1       = ca1.*contact1;
    sa1       = sa1.*contact1;
    contact2  = (l2<=sys.l0);
    ca2       = ca2.*contact2;
    sa2       = sa2.*contact2;
         
    fu = zeros(size(x, 1), size(u, 1), size(x, 2));
    fu(3,1,:) = ca1/sys.m;
    fu(3,2,:) = ca2/sys.m;
    fu(3,3,:) = sa1./(x(1,:)*sys.m);
    fu(3,4,:) = sa2./(l2*sys.m);
    
    fu(4,1,:) = sa1/sys.m;
    fu(4,2,:) = sa2/sys.m;
    fu(4,3,:) = -ca1./(x(1,:)*sys.m);
    fu(4,4,:) = -ca2./(l2*sys.m);
    
    fu(6,1,:) = sys.d*cos(x(2,:) - x(5,:));
    fu(6,2,:) = sys.d*cos(a2 - x(5,:));
    fu(6,3,:) = 1 + sys.d./x(1,:).*sin(x(2,:) - x(5,:));
    fu(6,4,:) = 1 + sys.d./l2.*sin(a2 - x(5,:));
end
