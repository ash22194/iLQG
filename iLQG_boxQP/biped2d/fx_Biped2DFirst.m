function fx = fx_Biped2DFirst(sys, x, u)
%FX_BIPED2D
%    DDX = FX_BIPED2D(sys, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    24-Apr-2020 16:17:49
l_point = sys.l_point;

U_DIMS_FREE = sys.U_DIMS_FREE;
U = zeros(4, size(u, 2));
U(U_DIMS_FREE, :) = u; 

X_DIMS_FREE = sys.X_DIMS_FREE;
X_DIMS_FIXED = sys.X_DIMS_FIXED;
X = zeros(6, size(x, 2));
X(X_DIMS_FREE, :) = x;
X(X_DIMS_FIXED, :) = repmat(l_point(X_DIMS_FIXED), [1, size(x,2)]);

x1 = X(1,:);
x2 = X(2,:);
x3 = X(3,:);
x4 = X(4,:);
x5 = X(5,:);
x6 = X(6,:);

u1 = U(1,:);
u2 = U(2,:);
u3 = U(3,:);
u4 = U(4,:);

m = sys.m;
l0 = sys.l0;
d = sys.d;
df = sys.df;
I = sys.I;
x_hip     = x1.*cos(x2);
z_hip     = x1.*sin(x2);
l2        = sqrt((x_hip + df).^2 + z_hip.^2);

u1 = u1.*(x1<=l0);
u3 = u3.*(x1<=l0);
u2 = u2.*(l2<=l0);
u4 = u4.*(l2<=l0);

t2 = cos(x2);
t3 = cos(x5);
t4 = sin(x2);
t5 = sin(x5);
t6 = x1.^2;
t8 = 1.0./I;
t10 = 1.0./m;
t11 = -x5;
t12 = 1.0./x1;
t7 = t4.^2;
t9 = t2.*x1;
t13 = 1.0./t6;
t14 = d.*t3.*x6;
t15 = d.*t5.*x6;
t17 = t11+x2;
t25 = t2.*t4.*t6.*2.0;
t16 = df+t9;
t18 = cos(t17);
t19 = sin(t17);
t20 = t7.*x1.*2.0;
t21 = t14+x3;
t22 = t15+x4;
t23 = t6.*t7;
t30 = -t25;
t24 = t16.^2;
t26 = t4.*t21;
t27 = t2.*t22;
t28 = d.*t19.*u1;
t29 = t2.*t16.*2.0;
t31 = t4.*t16.*x1.*2.0;
t32 = d.*t12.*t18.*u3;
t33 = t23+t24;
t34 = t20+t29;
t40 = t30+t31;
t35 = 1.0./t33;
t37 = sqrt(t33);
t36 = t35.^2;
t38 = 1.0./t37;
t44 = t24.*t35;
t47 = t29.*t35;
t49 = t31.*t35;
t50 = t2.*t16.*t35.*-2.0;
t51 = t4.*t16.*t35.*x1.*-2.0;
t39 = t38.^3;
t41 = t2.*t38;
t42 = t4.*t38.*x1;
t45 = t16.*t38;
t48 = -t44;
t61 = t24.*t34.*t36;
t64 = -t24.*t36.*(t25-t31);
t43 = -t42;
t46 = acos(t45);
t52 = t48+1.0;
t62 = (t16.*t34.*t39)./2.0;
t65 = t16.*t39.*(t25-t31).*(-1.0./2.0);
t68 = t50+t61;
t70 = t51+t64;
t53 = -t46;
t57 = sqrt(t52);
t63 = -t62;
t69 = t43+t65;
t54 = t53+x5;
t58 = 1.0./t57;
t67 = t41+t63;
t55 = cos(t54);
t56 = sin(t54);
t59 = d.*t56;
t60 = -t59;
t66 = t37+t60;

fx = zeros(size(X, 1), size(X, 1), size(X, 2));
fx(2,1,:) = t13.*(t26-t27);
fx(3,1,:) = -t10.*(-t41.*u2+t62.*u2+t4.*t13.*u3+(t34.*t39.*t57.*u4)./2.0-(t38.*t58.*u4.*(t61-t2.*t16.*t35.*2.0))./2.0);
fx(4,1,:) = t10.*(t2.*t13.*u3-t2.*t35.*u4+(t58.*u2.*(t61-t2.*t16.*t35.*2.0))./2.0+t16.*t34.*t36.*u4);
fx(6,1,:) = -t8.*(-t12.*u3-t38.*u4.*((t34.*t38)./2.0-d.*t55.*t58.*t67)+t13.*u3.*(x1+d.*t19)+(t34.*t39.*t66.*u4)./2.0+t58.*t59.*t67.*u2);
fx(1,2,:) = -t26+t27;
fx(2,2,:) = -t12.*(t2.*t21+t4.*t22);
fx(3,2,:) = -t10.*(t4.*u1+t42.*u2-t2.*t12.*u3+(t16.*t39.*u2.*(t25-t31))./2.0+(t39.*t57.*u4.*(t25-t31))./2.0-(t38.*t58.*u4.*(t49+t24.*t36.*(t25-t31)))./2.0);
fx(4,2,:) = t10.*(t2.*u1+(t58.*u2.*(t49+t24.*t36.*(t25-t31)))./2.0+t4.*t12.*u3+t4.*t35.*u4.*x1+t16.*t36.*u4.*(t25-t31));
fx(6,2,:) = t8.*(-t28+t32+t38.*u4.*((t38.*(t25-t31))./2.0+d.*t55.*t58.*(t42+(t16.*t39.*(t25-t31))./2.0))-(t39.*t66.*u4.*(t25-t31))./2.0+t58.*t59.*u2.*(t42+(t16.*t39.*(t25-t31))./2.0));
fx(1,3,:) = t2;
fx(2,3,:) = -t4.*t12;
fx(1,4,:) = t4;
fx(2,4,:) = t2.*t12;
fx(1,5,:) = -t2.*t15+t4.*t14;
fx(2,5,:) = t12.*(t2.*t14+t4.*t15);
fx(6,5,:) = -t8.*(-t28+t32+t59.*u2+d.*t38.*t55.*u4);
fx(1,6,:) = d.*t2.*t3+d.*t4.*t5;
fx(2,6,:) = t12.*(d.*t2.*t5-d.*t3.*t4);
fx(5,6,:) = ones(1,1,size(x,2));

fx = fx(X_DIMS_FREE, :, :);
fx = fx(:, X_DIMS_FREE, :);
end



