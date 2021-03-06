function out = dyn_gpu(sys, x, u)
%DYN
%    OUT = DYN(SYS, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    21-Sep-2020 10:35:56

m1 = sys.m(1);
m2 = sys.m(2);
m3 = sys.m(3);
l1 = sys.l(1);
l2 = sys.l(2);
l3 = sys.l(3);
g = sys.g;
X_DIMS_FREE = sys.X_DIMS_FREE; % Assuming they are sorted

th1 = x{1};
th2 = x{2};
th3 = x{3};
dth1 = x{4};
dth2 = x{5};
dth3 = x{6};

u1 = u{1};
u2 = u{2};
u3 = u{3};

out = cell(length(X_DIMS_FREE), 1);

for xxi=1:1:length(X_DIMS_FREE)
    xx = X_DIMS_FREE(xxi);
    if (xx==1)
        out{xxi} = dth1;
    elseif (xx==2)
        out{xxi} = dth2;
    elseif (xx==3)
        out{xxi} = dth3;
    elseif (xx==4)
        out{xxi} = arrayfun(@dynamics1, m1, m2, m3, l1, l2, l3, g, ...
                        th1, th2, th3, dth1, dth2, dth3, u1, u2, u3);
    elseif (xx==5)
        out{xxi} = arrayfun(@dynamics2, m1, m2, m3, l1, l2, l3, g, ...
                        th1, th2, th3, dth1, dth2, dth3, u1, u2, u3);
    elseif (xx==6)
        out{xxi} = arrayfun(@dynamics3, m1, m2, m3, l1, l2, l3, g, ...
                        th1, th2, th3, dth1, dth2, dth3, u1, u2, u3);
    end
end

    function d_dth1 = dynamics1(m1, m2, m3, l1, l2, l3, g, ...
                            th1, th2, th3, dth1, dth2, dth3, ...
                            u1, u2, u3)

        t2 = cos(th2);
        t3 = cos(th3);
        t4 = sin(th1);
        t5 = sin(th2);
        t6 = sin(th3);
        t7 = th1+th2;
        t8 = th2+th3;
        t9 = dth1.^2;
        t10 = dth2.^2;
        t11 = dth3.^2;
        t12 = l1.^2;
        t13 = l1.^3;
        t14 = l2.^2;
        t15 = l2.^3;
        t16 = l3.^2;
        t17 = l3.^3;
        t18 = m2.^2;
        t19 = m3.^2;
        t20 = m3.^3;
        t21 = th2.*2.0;
        t22 = th3.*2.0;
        t33 = 1.0./l3;
        t34 = -th1;
        t35 = -th2;
        t36 = -th3;
        t39 = m1.*m2.*1.6e+1;
        t40 = m1.*m3.*3.0e+1;
        t54 = m2.*m3.*7.5e+1;
        t23 = cos(t21);
        t24 = cos(t22);
        t25 = sin(t21);
        t26 = sin(t22);
        t27 = cos(t8);
        t28 = sin(t7);
        t29 = sin(t8);
        t30 = t7+th3;
        t31 = 1.0./t12;
        t32 = 1.0./t14;
        t37 = -t22;
        t41 = t7+th2;
        t42 = t22+th1;
        t43 = t8+th3;
        t44 = t8+th2;
        t51 = t7+t22;
        t52 = t19.*1.8e+1;
        t53 = t18.*3.0e+1;
        t55 = t35+th1;
        t57 = t36+th2;
        t58 = t8.*2.0;
        t63 = t7+t36;
        t65 = t8+t34;
        t38 = sin(t30);
        t45 = cos(t43);
        t46 = cos(t44);
        t47 = sin(t41);
        t48 = sin(t42);
        t49 = sin(t43);
        t50 = sin(t44);
        t56 = t37+th1;
        t59 = cos(t57);
        t60 = sin(t55);
        t62 = sin(t57);
        t64 = t55+th3;
        t66 = sin(t51);
        t67 = cos(t58);
        t68 = sin(t58);
        t69 = t8+t30;
        t70 = sin(t63);
        t72 = sin(t65);
        t73 = m1.*m3.*t24.*1.8e+1;
        t74 = m2.*m3.*t24.*2.7e+1;
        t75 = m2.*m3.*t23.*4.5e+1;
        t76 = t34+t43;
        t78 = t18.*t23.*1.8e+1;
        t79 = t23.*t52;
        t85 = t19.*t23.*-1.8e+1;
        t61 = sin(t56);
        t71 = sin(t64);
        t77 = sin(t69);
        t80 = -t73;
        t81 = -t74;
        t82 = -t75;
        t83 = sin(t76);
        t84 = -t78;
        t86 = m2.*m3.*t67.*9.0;
        t87 = t39+t40+t52+t53+t54+t80+t81+t82+t84+t85+t86;
        t88 = 1.0./t87;
        
        d_dth1 = (t31.*t33.*t88.*(l2.*l3.*m2.*u1.*8.0-l2.*l3.*m2.*u2.*8.0+l2.*l3.*m3.*u1.*1.5e+1-l2.*l3.*m3.*u2.*1.5e+1-l1.*l3.*m2.*t2.*u2.*1.2e+1+l1.*l3.*m2.*t2.*u3.*1.2e+1-l1.*l3.*m3.*t2.*u2.*1.5e+1+l1.*l3.*m3.*t2.*u3.*1.5e+1-l2.*l3.*m3.*t24.*u1.*9.0+l2.*l3.*m3.*t24.*u2.*9.0-l1.*l2.*m2.*t27.*u3.*3.0-l1.*l2.*m3.*t27.*u3.*1.8e+1+l1.*l3.*m3.*t45.*u2.*9.0-l1.*l3.*m3.*t45.*u3.*9.0+l1.*l2.*m2.*t59.*u3.*9.0+l1.*l2.*m3.*t59.*u3.*1.8e+1+l1.*l3.*t5.*t9.*t14.*t18.*4.0+l1.*l3.*t5.*t9.*t14.*t19.*6.0+l1.*l3.*t5.*t10.*t14.*t18.*4.0+l1.*l3.*t5.*t10.*t14.*t19.*6.0+l2.*l3.*t9.*t12.*t18.*t25.*3.0+l2.*l3.*t9.*t12.*t19.*t25.*3.0+l1.*l2.*t9.*t16.*t19.*t29.*(3.0./2.0)+l1.*l2.*t10.*t16.*t19.*t29.*(3.0./2.0)+l1.*l2.*t11.*t16.*t19.*t29.*(3.0./2.0)+l1.*l2.*t9.*t16.*t19.*t62.*(3.0./2.0)+l1.*l2.*t10.*t16.*t19.*t62.*(3.0./2.0)+l1.*l2.*t11.*t16.*t19.*t62.*(3.0./2.0)-g.*l1.*l2.*l3.*t4.*t18.*5.0-g.*l1.*l2.*l3.*t4.*t19.*3.0+g.*l1.*l2.*l3.*t18.*t47.*3.0+g.*l1.*l2.*l3.*t19.*t47.*3.0-g.*l1.*l2.*l3.*m1.*m2.*t4.*4.0-g.*l1.*l2.*l3.*m1.*m3.*t4.*(1.5e+1./2.0)-g.*l1.*l2.*l3.*m2.*m3.*t4.*(2.5e+1./2.0)+g.*l1.*l2.*l3.*m1.*m3.*t48.*(9.0./4.0)+g.*l1.*l2.*l3.*m2.*m3.*t47.*(1.5e+1./2.0)+g.*l1.*l2.*l3.*m2.*m3.*t48.*(9.0./4.0)+g.*l1.*l2.*l3.*m1.*m3.*t61.*(9.0./4.0)+g.*l1.*l2.*l3.*m2.*m3.*t61.*(9.0./4.0)-g.*l1.*l2.*l3.*m2.*m3.*t77.*(3.0./2.0)+dth1.*dth2.*l1.*l3.*t5.*t14.*t18.*8.0+dth1.*dth2.*l1.*l3.*t5.*t14.*t19.*1.2e+1+dth1.*dth2.*l1.*l2.*t16.*t19.*t29.*3.0+dth1.*dth3.*l1.*l2.*t16.*t19.*t29.*3.0+dth2.*dth3.*l1.*l2.*t16.*t19.*t29.*3.0+dth1.*dth2.*l1.*l2.*t16.*t19.*t62.*3.0+dth1.*dth3.*l1.*l2.*t16.*t19.*t62.*3.0+dth2.*dth3.*l1.*l2.*t16.*t19.*t62.*3.0+l1.*l3.*m2.*m3.*t5.*t9.*t14.*(2.5e+1./2.0)+l1.*l3.*m2.*m3.*t5.*t10.*t14.*(2.5e+1./2.0)+l2.*l3.*m2.*m3.*t9.*t12.*t25.*(1.5e+1./2.0)+l1.*l2.*m2.*m3.*t9.*t16.*t29+l1.*l2.*m2.*m3.*t10.*t16.*t29+l1.*l2.*m2.*m3.*t11.*t16.*t29-l1.*l3.*m2.*m3.*t9.*t14.*t49.*(3.0./2.0)-l1.*l3.*m2.*m3.*t10.*t14.*t49.*(3.0./2.0)+l1.*l2.*m2.*m3.*t9.*t16.*t62.*3.0+l1.*l2.*m2.*m3.*t10.*t16.*t62.*3.0+l1.*l2.*m2.*m3.*t11.*t16.*t62.*3.0-l2.*l3.*m2.*m3.*t9.*t12.*t68.*(3.0./2.0)+dth1.*dth2.*l1.*l3.*m2.*m3.*t5.*t14.*2.5e+1+dth1.*dth2.*l1.*l2.*m2.*m3.*t16.*t29.*2.0+dth1.*dth3.*l1.*l2.*m2.*m3.*t16.*t29.*2.0+dth2.*dth3.*l1.*l2.*m2.*m3.*t16.*t29.*2.0-dth1.*dth2.*l1.*l3.*m2.*m3.*t14.*t49.*3.0+dth1.*dth2.*l1.*l2.*m2.*m3.*t16.*t62.*6.0+dth1.*dth3.*l1.*l2.*m2.*m3.*t16.*t62.*6.0+dth2.*dth3.*l1.*l2.*m2.*m3.*t16.*t62.*6.0).*6.0)./l2;

    end

    function d_dth2 = dynamics2(m1, m2, m3, l1, l2, l3, g, ...
                            th1, th2, th3, dth1, dth2, dth3, ...
                            u1, u2, u3)
        t2 = cos(th2);
        t3 = cos(th3);
        t4 = sin(th1);
        t5 = sin(th2);
        t6 = sin(th3);
        t7 = th1+th2;
        t8 = th2+th3;
        t9 = dth1.^2;
        t10 = dth2.^2;
        t11 = dth3.^2;
        t12 = l1.^2;
        t13 = l1.^3;
        t14 = l2.^2;
        t15 = l2.^3;
        t16 = l3.^2;
        t17 = l3.^3;
        t18 = m2.^2;
        t19 = m3.^2;
        t20 = m3.^3;
        t21 = th2.*2.0;
        t22 = th3.*2.0;
        t33 = 1.0./l3;
        t34 = -th1;
        t35 = -th2;
        t36 = -th3;
        t39 = m1.*m2.*1.6e+1;
        t40 = m1.*m3.*3.0e+1;
        t54 = m2.*m3.*7.5e+1;
        t23 = cos(t21);
        t24 = cos(t22);
        t25 = sin(t21);
        t26 = sin(t22);
        t27 = cos(t8);
        t28 = sin(t7);
        t29 = sin(t8);
        t30 = t7+th3;
        t31 = 1.0./t12;
        t32 = 1.0./t14;
        t37 = -t22;
        t41 = t7+th2;
        t42 = t22+th1;
        t43 = t8+th3;
        t44 = t8+th2;
        t51 = t7+t22;
        t52 = t19.*1.8e+1;
        t53 = t18.*3.0e+1;
        t55 = t35+th1;
        t57 = t36+th2;
        t58 = t8.*2.0;
        t63 = t7+t36;
        t65 = t8+t34;
        t38 = sin(t30);
        t45 = cos(t43);
        t46 = cos(t44);
        t47 = sin(t41);
        t48 = sin(t42);
        t49 = sin(t43);
        t50 = sin(t44);
        t56 = t37+th1;
        t59 = cos(t57);
        t60 = sin(t55);
        t62 = sin(t57);
        t64 = t55+th3;
        t66 = sin(t51);
        t67 = cos(t58);
        t68 = sin(t58);
        t69 = t8+t30;
        t70 = sin(t63);
        t72 = sin(t65);
        t73 = m1.*m3.*t24.*1.8e+1;
        t74 = m2.*m3.*t24.*2.7e+1;
        t75 = m2.*m3.*t23.*4.5e+1;
        t76 = t34+t43;
        t78 = t18.*t23.*1.8e+1;
        t79 = t23.*t52;
        t85 = t19.*t23.*-1.8e+1;
        t61 = sin(t56);
        t71 = sin(t64);
        t77 = sin(t69);
        t80 = -t73;
        t81 = -t74;
        t82 = -t75;
        t83 = sin(t76);
        t84 = -t78;
        t86 = m2.*m3.*t67.*9.0;
        t87 = t39+t40+t52+t53+t54+t80+t81+t82+t84+t85+t86;
        t88 = 1.0./t87;
        
        d_dth2 = t31.*t32.*t33.*t88.*(l3.*m1.*t12.*u2.*-8.0+l3.*m1.*t12.*u3.*8.0-l3.*m2.*t12.*u2.*2.4e+1+l3.*m2.*t12.*u3.*2.4e+1+l3.*m2.*t14.*u1.*8.0-l3.*m3.*t12.*u2.*1.5e+1-l3.*m2.*t14.*u2.*8.0+l3.*m3.*t12.*u3.*1.5e+1+l3.*m3.*t14.*u1.*1.5e+1-l3.*m3.*t14.*u2.*1.5e+1+l2.*m1.*t3.*t12.*u3.*1.2e+1+l2.*m2.*t3.*t12.*u3.*2.7e+1+l2.*m3.*t3.*t12.*u3.*1.8e+1-l3.*m3.*t14.*t24.*u1.*9.0+l3.*m3.*t14.*t24.*u2.*9.0-l1.*m2.*t14.*t27.*u3.*3.0-l1.*m3.*t14.*t27.*u3.*1.8e+1-l2.*m2.*t12.*t46.*u3.*9.0-l2.*m3.*t12.*t46.*u3.*1.8e+1+l1.*m2.*t14.*t59.*u3.*9.0+l1.*m3.*t14.*t59.*u3.*1.8e+1+l3.*m3.*t12.*t67.*u2.*9.0-l3.*m3.*t12.*t67.*u3.*9.0+l2.*l3.*t5.*t9.*t13.*t18.*1.2e+1+l1.*l3.*t5.*t9.*t15.*t18.*4.0+l2.*l3.*t5.*t9.*t13.*t19.*6.0+l1.*l3.*t5.*t9.*t15.*t19.*6.0+l1.*l3.*t5.*t10.*t15.*t18.*4.0+l1.*l3.*t5.*t10.*t15.*t19.*6.0-l2.*t6.*t9.*t12.*t16.*t19.*(3.0./2.0)-l2.*t6.*t10.*t12.*t16.*t19.*(3.0./2.0)-l2.*t6.*t11.*t12.*t16.*t19.*(3.0./2.0)+l3.*t9.*t12.*t14.*t18.*t25.*6.0+l3.*t9.*t12.*t14.*t19.*t25.*6.0+l3.*t10.*t12.*t14.*t18.*t25.*3.0+l3.*t10.*t12.*t14.*t19.*t25.*3.0+l1.*t9.*t14.*t16.*t19.*t29.*(3.0./2.0)+l1.*t10.*t14.*t16.*t19.*t29.*(3.0./2.0)+l1.*t11.*t14.*t16.*t19.*t29.*(3.0./2.0)+l2.*t9.*t12.*t16.*t19.*t50.*(3.0./2.0)+l2.*t10.*t12.*t16.*t19.*t50.*(3.0./2.0)+l2.*t11.*t12.*t16.*t19.*t50.*(3.0./2.0)+l1.*t9.*t14.*t16.*t19.*t62.*(3.0./2.0)+l1.*t10.*t14.*t16.*t19.*t62.*(3.0./2.0)+l1.*t11.*t14.*t16.*t19.*t62.*(3.0./2.0)+l1.*l2.*l3.*m2.*t2.*u1.*1.2e+1-l1.*l2.*l3.*m2.*t2.*u2.*2.4e+1+l1.*l2.*l3.*m3.*t2.*u1.*1.5e+1+l1.*l2.*l3.*m2.*t2.*u3.*1.2e+1-l1.*l2.*l3.*m3.*t2.*u2.*3.0e+1+l1.*l2.*l3.*m3.*t2.*u3.*1.5e+1-l1.*l2.*l3.*m3.*t45.*u1.*9.0+l1.*l2.*l3.*m3.*t45.*u2.*1.8e+1-l1.*l2.*l3.*m3.*t45.*u3.*9.0-g.*l1.*l3.*t4.*t14.*t18.*5.0-g.*l1.*l3.*t4.*t14.*t19.*3.0+g.*l2.*l3.*t12.*t18.*t28.*6.0+g.*l2.*l3.*t12.*t19.*t28.*3.0+g.*l1.*l3.*t14.*t18.*t47.*3.0+g.*l1.*l3.*t14.*t19.*t47.*3.0-g.*l2.*l3.*t12.*t18.*t60.*6.0-g.*l2.*l3.*t12.*t19.*t60.*3.0+dth1.*dth2.*l1.*l3.*t5.*t15.*t18.*8.0+dth1.*dth2.*l1.*l3.*t5.*t15.*t19.*1.2e+1-g.*l1.*l3.*m1.*m2.*t4.*t14.*4.0-g.*l1.*l3.*m1.*m3.*t4.*t14.*(1.5e+1./2.0)-g.*l1.*l3.*m2.*m3.*t4.*t14.*(2.5e+1./2.0)+g.*l2.*l3.*m1.*m2.*t12.*t28+g.*l2.*l3.*m1.*m3.*t12.*t28.*(5.0./4.0)+g.*l2.*l3.*m2.*m3.*t12.*t28.*(4.5e+1./4.0)+g.*l1.*l3.*m1.*m3.*t14.*t48.*(9.0./4.0)+g.*l1.*l3.*m2.*m3.*t14.*t47.*(1.5e+1./2.0)+g.*l1.*l3.*m2.*m3.*t14.*t48.*(9.0./4.0)-g.*l2.*l3.*m1.*m2.*t12.*t60.*3.0-g.*l2.*l3.*m1.*m3.*t12.*t60.*(1.5e+1./4.0)-g.*l2.*l3.*m2.*m3.*t12.*t60.*(4.5e+1./4.0)+g.*l1.*l3.*m1.*m3.*t14.*t61.*(9.0./4.0)+g.*l1.*l3.*m2.*m3.*t14.*t61.*(9.0./4.0)-g.*l2.*l3.*m1.*m3.*t12.*t66.*(3.0./4.0)-g.*l2.*l3.*m2.*m3.*t12.*t66.*(9.0./4.0)-g.*l1.*l3.*m2.*m3.*t14.*t77.*(3.0./2.0)-g.*l2.*l3.*m1.*m3.*t12.*t83.*(9.0./4.0)-g.*l2.*l3.*m2.*m3.*t12.*t83.*(9.0./4.0)-dth1.*dth2.*l2.*t6.*t12.*t16.*t19.*3.0-dth1.*dth3.*l2.*t6.*t12.*t16.*t19.*3.0-dth2.*dth3.*l2.*t6.*t12.*t16.*t19.*3.0+dth1.*dth2.*l3.*t12.*t14.*t18.*t25.*6.0+dth1.*dth2.*l3.*t12.*t14.*t19.*t25.*6.0+dth1.*dth2.*l1.*t14.*t16.*t19.*t29.*3.0+dth1.*dth3.*l1.*t14.*t16.*t19.*t29.*3.0+dth2.*dth3.*l1.*t14.*t16.*t19.*t29.*3.0+dth1.*dth2.*l2.*t12.*t16.*t19.*t50.*3.0+dth1.*dth3.*l2.*t12.*t16.*t19.*t50.*3.0+dth2.*dth3.*l2.*t12.*t16.*t19.*t50.*3.0+dth1.*dth2.*l1.*t14.*t16.*t19.*t62.*3.0+dth1.*dth3.*l1.*t14.*t16.*t19.*t62.*3.0+dth2.*dth3.*l1.*t14.*t16.*t19.*t62.*3.0+l2.*l3.*m1.*m2.*t5.*t9.*t13.*4.0+l2.*l3.*m1.*m3.*t5.*t9.*t13.*5.0+l2.*l3.*m2.*m3.*t5.*t9.*t13.*(4.5e+1./2.0)+l1.*l3.*m2.*m3.*t5.*t9.*t15.*(2.5e+1./2.0)+l1.*l3.*m2.*m3.*t5.*t10.*t15.*(2.5e+1./2.0)-l2.*l3.*m1.*m3.*t9.*t13.*t49.*3.0-l2.*l3.*m2.*m3.*t9.*t13.*t49.*(9.0./2.0)-l1.*l3.*m2.*m3.*t9.*t15.*t49.*(3.0./2.0)-l1.*l3.*m2.*m3.*t10.*t15.*t49.*(3.0./2.0)-l2.*m1.*m3.*t6.*t9.*t12.*t16.*4.0-l2.*m1.*m3.*t6.*t10.*t12.*t16.*4.0-l2.*m2.*m3.*t6.*t9.*t12.*t16.*9.0-l2.*m1.*m3.*t6.*t11.*t12.*t16.*4.0-l2.*m2.*m3.*t6.*t10.*t12.*t16.*9.0-l2.*m2.*m3.*t6.*t11.*t12.*t16.*9.0-l3.*m1.*m3.*t9.*t12.*t14.*t26.*3.0+l3.*m2.*m3.*t9.*t12.*t14.*t25.*1.5e+1-l3.*m1.*m3.*t10.*t12.*t14.*t26.*3.0-l3.*m2.*m3.*t9.*t12.*t14.*t26.*(9.0./2.0)+l3.*m2.*m3.*t10.*t12.*t14.*t25.*(1.5e+1./2.0)-l3.*m2.*m3.*t10.*t12.*t14.*t26.*(9.0./2.0)+l1.*m2.*m3.*t9.*t14.*t16.*t29+l1.*m2.*m3.*t10.*t14.*t16.*t29+l1.*m2.*m3.*t11.*t14.*t16.*t29+l2.*m2.*m3.*t9.*t12.*t16.*t50.*3.0+l2.*m2.*m3.*t10.*t12.*t16.*t50.*3.0+l2.*m2.*m3.*t11.*t12.*t16.*t50.*3.0+l1.*m2.*m3.*t9.*t14.*t16.*t62.*3.0+l1.*m2.*m3.*t10.*t14.*t16.*t62.*3.0+l1.*m2.*m3.*t11.*t14.*t16.*t62.*3.0-l3.*m2.*m3.*t9.*t12.*t14.*t68.*(3.0./2.0)+dth1.*dth2.*l1.*l3.*m2.*m3.*t5.*t15.*2.5e+1-dth1.*dth2.*l1.*l3.*m2.*m3.*t15.*t49.*3.0-dth1.*dth2.*l2.*m1.*m3.*t6.*t12.*t16.*8.0-dth1.*dth2.*l2.*m2.*m3.*t6.*t12.*t16.*1.8e+1-dth1.*dth3.*l2.*m1.*m3.*t6.*t12.*t16.*8.0-dth1.*dth3.*l2.*m2.*m3.*t6.*t12.*t16.*1.8e+1-dth2.*dth3.*l2.*m1.*m3.*t6.*t12.*t16.*8.0-dth2.*dth3.*l2.*m2.*m3.*t6.*t12.*t16.*1.8e+1-dth1.*dth2.*l3.*m1.*m3.*t12.*t14.*t26.*6.0+dth1.*dth2.*l3.*m2.*m3.*t12.*t14.*t25.*1.5e+1-dth1.*dth2.*l3.*m2.*m3.*t12.*t14.*t26.*9.0+dth1.*dth2.*l1.*m2.*m3.*t14.*t16.*t29.*2.0+dth1.*dth3.*l1.*m2.*m3.*t14.*t16.*t29.*2.0+dth2.*dth3.*l1.*m2.*m3.*t14.*t16.*t29.*2.0+dth1.*dth2.*l2.*m2.*m3.*t12.*t16.*t50.*6.0+dth1.*dth3.*l2.*m2.*m3.*t12.*t16.*t50.*6.0+dth2.*dth3.*l2.*m2.*m3.*t12.*t16.*t50.*6.0+dth1.*dth2.*l1.*m2.*m3.*t14.*t16.*t62.*6.0+dth1.*dth3.*l1.*m2.*m3.*t14.*t16.*t62.*6.0+dth2.*dth3.*l1.*m2.*m3.*t14.*t16.*t62.*6.0).*-6.0;

    end

    function d_dth3 = dynamics3(m1, m2, m3, l1, l2, l3, g, ...
                                th1, th2, th3, dth1, dth2, dth3, ...
                                u1, u2, u3)
        t2 = cos(th2);
        t3 = cos(th3);
        t4 = sin(th1);
        t5 = sin(th2);
        t6 = sin(th3);
        t7 = th1+th2;
        t8 = th2+th3;
        t9 = dth1.^2;
        t10 = dth2.^2;
        t11 = dth3.^2;
        t12 = l1.^2;
        t13 = l1.^3;
        t14 = l2.^2;
        t15 = l2.^3;
        t16 = l3.^2;
        t17 = l3.^3;
        t18 = m2.^2;
        t19 = m3.^2;
        t20 = m3.^3;
        t21 = th2.*2.0;
        t22 = th3.*2.0;
        t33 = 1.0./l3;
        t34 = -th1;
        t35 = -th2;
        t36 = -th3;
        t39 = m1.*m2.*1.6e+1;
        t40 = m1.*m3.*3.0e+1;
        t54 = m2.*m3.*7.5e+1;
        t23 = cos(t21);
        t24 = cos(t22);
        t25 = sin(t21);
        t26 = sin(t22);
        t27 = cos(t8);
        t28 = sin(t7);
        t29 = sin(t8);
        t30 = t7+th3;
        t31 = 1.0./t12;
        t32 = 1.0./t14;
        t37 = -t22;
        t41 = t7+th2;
        t42 = t22+th1;
        t43 = t8+th3;
        t44 = t8+th2;
        t51 = t7+t22;
        t52 = t19.*1.8e+1;
        t53 = t18.*3.0e+1;
        t55 = t35+th1;
        t57 = t36+th2;
        t58 = t8.*2.0;
        t63 = t7+t36;
        t65 = t8+t34;
        t38 = sin(t30);
        t45 = cos(t43);
        t46 = cos(t44);
        t47 = sin(t41);
        t48 = sin(t42);
        t49 = sin(t43);
        t50 = sin(t44);
        t56 = t37+th1;
        t59 = cos(t57);
        t60 = sin(t55);
        t62 = sin(t57);
        t64 = t55+th3;
        t66 = sin(t51);
        t67 = cos(t58);
        t68 = sin(t58);
        t69 = t8+t30;
        t70 = sin(t63);
        t72 = sin(t65);
        t73 = m1.*m3.*t24.*1.8e+1;
        t74 = m2.*m3.*t24.*2.7e+1;
        t75 = m2.*m3.*t23.*4.5e+1;
        t76 = t34+t43;
        t78 = t18.*t23.*1.8e+1;
        t79 = t23.*t52;
        t85 = t19.*t23.*-1.8e+1;
        t61 = sin(t56);
        t71 = sin(t64);
        t77 = sin(t69);
        t80 = -t73;
        t81 = -t74;
        t82 = -t75;
        t83 = sin(t76);
        t84 = -t78;
        t86 = m2.*m3.*t67.*9.0;
        t87 = t39+t40+t52+t53+t54+t80+t81+t82+t84+t85+t86;
        t88 = 1.0./t87;
        
        d_dth3 = (t32.*t88.*(l1.*t14.*t18.*u3.*-1.5e+1-l1.*t14.*t19.*u3.*3.6e+1+l1.*t16.*t19.*u2.*1.5e+1-l1.*t16.*t19.*u3.*1.5e+1-l1.*m1.*m2.*t14.*u3.*8.0-l1.*m1.*m3.*t14.*u3.*2.4e+1+l1.*m1.*m3.*t16.*u2.*8.0-l1.*m2.*m3.*t14.*u3.*6.0e+1-l1.*m1.*m3.*t16.*u3.*8.0+l1.*m2.*m3.*t16.*u2.*2.4e+1-l1.*m2.*m3.*t16.*u3.*2.4e+1-l2.*t2.*t16.*t19.*u1.*1.5e+1+l2.*t2.*t16.*t19.*u2.*1.5e+1+l1.*t14.*t18.*t23.*u3.*9.0+l1.*t14.*t19.*t23.*u3.*3.6e+1-l3.*t14.*t19.*t27.*u2.*1.8e+1+l2.*t16.*t19.*t45.*u1.*9.0-l2.*t16.*t19.*t45.*u2.*9.0-l3.*t14.*t19.*t59.*u1.*1.8e+1+l3.*t14.*t27.*t52.*u1-l1.*t16.*t19.*t67.*u2.*9.0+l1.*t16.*t19.*t67.*u3.*9.0+l3.*t14.*t52.*t59.*u2-l1.*l2.*l3.*t3.*t19.*u3.*3.6e+1+l1.*l2.*l3.*t3.*t52.*u2-l1.*l2.*l3.*t19.*t46.*u2.*1.8e+1+l1.*l2.*l3.*t19.*t46.*u3.*3.6e+1-l2.*m2.*m3.*t2.*t16.*u1.*1.2e+1+l2.*m2.*m3.*t2.*t16.*u2.*1.2e+1+l1.*m2.*m3.*t14.*t23.*u3.*3.6e+1+l3.*m2.*m3.*t14.*t27.*u1.*3.0-l3.*m2.*m3.*t14.*t27.*u2.*3.0-l3.*m2.*m3.*t14.*t59.*u1.*9.0+l3.*m2.*m3.*t14.*t59.*u2.*9.0+l1.*l2.*t6.*t9.*t17.*t20.*(3.0./2.0)+l1.*l2.*t6.*t10.*t17.*t20.*(3.0./2.0)+l1.*l2.*t6.*t11.*t17.*t20.*(3.0./2.0)-l1.*l2.*t9.*t17.*t20.*t50.*(3.0./2.0)-l1.*l2.*t10.*t17.*t20.*t50.*(3.0./2.0)-l1.*l2.*t11.*t17.*t20.*t50.*(3.0./2.0)-l2.*t5.*t9.*t12.*t16.*t20.*6.0-l1.*t9.*t14.*t16.*t20.*t25.*3.0-l1.*t10.*t14.*t16.*t20.*t25.*3.0-g.*l1.*l2.*t16.*t20.*t28.*3.0+g.*l1.*l2.*t16.*t20.*t60.*3.0+dth1.*dth2.*l1.*l2.*t6.*t17.*t20.*3.0+dth1.*dth3.*l1.*l2.*t6.*t17.*t20.*3.0+dth2.*dth3.*l1.*l2.*t6.*t17.*t20.*3.0-dth1.*dth2.*l1.*l2.*t17.*t20.*t50.*3.0-dth1.*dth3.*l1.*l2.*t17.*t20.*t50.*3.0-dth2.*dth3.*l1.*l2.*t17.*t20.*t50.*3.0-dth1.*dth2.*l1.*t14.*t16.*t20.*t25.*6.0+l1.*l2.*l3.*m1.*m3.*t3.*u2.*1.2e+1-l1.*l2.*l3.*m1.*m3.*t3.*u3.*2.4e+1+l1.*l2.*l3.*m2.*m3.*t3.*u2.*2.7e+1-l1.*l2.*l3.*m2.*m3.*t3.*u3.*5.4e+1-l1.*l2.*l3.*m2.*m3.*t46.*u2.*9.0+l1.*l2.*l3.*m2.*m3.*t46.*u3.*1.8e+1-g.*l1.*l2.*m1.*t16.*t19.*t28.*(5.0./4.0)-g.*l1.*l2.*m2.*t16.*t19.*t28.*(4.5e+1./4.0)-g.*l1.*l2.*m3.*t16.*t18.*t28.*6.0+g.*l1.*l3.*m1.*t14.*t19.*t38.*(3.0./2.0)+g.*l1.*l3.*m2.*t14.*t19.*t38.*(3.0./2.0)-g.*l1.*l3.*m3.*t14.*t18.*t38.*(3.0./4.0)+g.*l1.*l2.*m1.*t16.*t19.*t60.*(1.5e+1./4.0)+g.*l1.*l2.*m2.*t16.*t19.*t60.*(4.5e+1./4.0)+g.*l1.*l2.*m3.*t16.*t18.*t60.*6.0+g.*l1.*l2.*m1.*t16.*t19.*t66.*(3.0./4.0)+g.*l1.*l2.*m2.*t16.*t19.*t66.*(9.0./4.0)-g.*l1.*l3.*m1.*t14.*t19.*t70.*(3.0./2.0)+g.*l1.*l3.*m1.*t14.*t19.*t71.*(9.0./2.0)-g.*l1.*l3.*m2.*t14.*t19.*t70.*(9.0./2.0)-g.*l1.*l3.*m3.*t14.*t18.*t70.*(9.0./4.0)+g.*l1.*l3.*m1.*t14.*t19.*t72.*(9.0./2.0)+g.*l1.*l3.*m2.*t14.*t19.*t71.*(9.0./2.0)+g.*l1.*l3.*m3.*t14.*t18.*t71.*(9.0./4.0)+g.*l1.*l3.*m2.*t14.*t19.*t72.*(3.0./2.0)-g.*l1.*l3.*m3.*t14.*t18.*t72.*(3.0./4.0)+g.*l1.*l2.*m1.*t16.*t19.*t83.*(9.0./4.0)+g.*l1.*l2.*m2.*t16.*t19.*t83.*(9.0./4.0)+l1.*l3.*m1.*t6.*t9.*t15.*t19.*1.2e+1+l1.*l2.*m1.*t6.*t9.*t17.*t19.*4.0+l1.*l3.*m1.*t6.*t10.*t15.*t19.*1.2e+1+l1.*l3.*m2.*t6.*t9.*t15.*t19.*1.5e+1+l1.*l3.*m3.*t6.*t9.*t15.*t18.*(9.0./2.0)+l1.*l2.*m1.*t6.*t10.*t17.*t19.*4.0+l1.*l2.*m2.*t6.*t9.*t17.*t19.*9.0+l1.*l3.*m2.*t6.*t10.*t15.*t19.*1.5e+1+l1.*l3.*m3.*t6.*t10.*t15.*t18.*(9.0./2.0)+l1.*l2.*m1.*t6.*t11.*t17.*t19.*4.0+l1.*l2.*m2.*t6.*t10.*t17.*t19.*9.0+l1.*l2.*m2.*t6.*t11.*t17.*t19.*9.0-l1.*l3.*m2.*t9.*t15.*t19.*t50.*3.0-l1.*l3.*m3.*t9.*t15.*t18.*t50.*(3.0./2.0)-l1.*l2.*m2.*t9.*t17.*t19.*t50.*3.0-l1.*l3.*m2.*t10.*t15.*t19.*t50.*3.0-l1.*l3.*m3.*t10.*t15.*t18.*t50.*(3.0./2.0)-l1.*l2.*m2.*t10.*t17.*t19.*t50.*3.0-l1.*l2.*m2.*t11.*t17.*t19.*t50.*3.0-l2.*m1.*t5.*t9.*t12.*t16.*t19.*5.0-l2.*m2.*t5.*t9.*t12.*t16.*t19.*(4.5e+1./2.0)-l2.*m3.*t5.*t9.*t12.*t16.*t18.*1.2e+1+l1.*m1.*t9.*t14.*t16.*t19.*t26.*6.0-l1.*m2.*t9.*t14.*t16.*t19.*t25.*(1.5e+1./2.0)-l1.*m3.*t9.*t14.*t16.*t18.*t25.*3.0+l1.*m1.*t10.*t14.*t16.*t19.*t26.*6.0+l1.*m2.*t9.*t14.*t16.*t19.*t26.*9.0-l1.*m2.*t10.*t14.*t16.*t19.*t25.*(1.5e+1./2.0)-l1.*m3.*t10.*t14.*t16.*t18.*t25.*3.0+l3.*m1.*t9.*t12.*t14.*t19.*t29.*6.0+l1.*m1.*t11.*t14.*t16.*t19.*t26.*3.0+l1.*m2.*t10.*t14.*t16.*t19.*t26.*9.0+l3.*m2.*t9.*t12.*t14.*t19.*t29.*3.0-l3.*m3.*t9.*t12.*t14.*t18.*t29.*(3.0./2.0)+l1.*m2.*t11.*t14.*t16.*t19.*t26.*(9.0./2.0)+l2.*m1.*t9.*t12.*t16.*t19.*t49.*3.0+l2.*m2.*t9.*t12.*t16.*t19.*t49.*(9.0./2.0)-l3.*m1.*t9.*t12.*t14.*t19.*t62.*6.0-l3.*m2.*t9.*t12.*t14.*t19.*t62.*9.0-l3.*m3.*t9.*t12.*t14.*t18.*t62.*(9.0./2.0)-l1.*m2.*t9.*t14.*t16.*t19.*t68.*(3.0./2.0)-l1.*m2.*t10.*t14.*t16.*t19.*t68.*(3.0./2.0)-l1.*m2.*t11.*t14.*t16.*t19.*t68.*(3.0./2.0)+dth1.*dth2.*l1.*l3.*m1.*t6.*t15.*t19.*2.4e+1+dth1.*dth2.*l1.*l2.*m1.*t6.*t17.*t19.*8.0+dth1.*dth2.*l1.*l3.*m2.*t6.*t15.*t19.*3.0e+1+dth1.*dth2.*l1.*l3.*m3.*t6.*t15.*t18.*9.0+dth1.*dth3.*l1.*l2.*m1.*t6.*t17.*t19.*8.0+dth2.*dth3.*l1.*l2.*m1.*t6.*t17.*t19.*8.0+dth1.*dth2.*l1.*l2.*m2.*t6.*t17.*t52+dth1.*dth3.*l1.*l2.*m2.*t6.*t17.*t52+dth2.*dth3.*l1.*l2.*m2.*t6.*t17.*t52-dth1.*dth2.*l1.*l3.*m2.*t15.*t19.*t50.*6.0-dth1.*dth2.*l1.*l3.*m3.*t15.*t18.*t50.*3.0-dth1.*dth2.*l1.*l2.*m2.*t17.*t19.*t50.*6.0-dth1.*dth3.*l1.*l2.*m2.*t17.*t19.*t50.*6.0-dth2.*dth3.*l1.*l2.*m2.*t17.*t19.*t50.*6.0-g.*l1.*l2.*m1.*m2.*m3.*t16.*t28+(g.*l1.*l3.*m1.*m2.*m3.*t14.*t38)./4.0+g.*l1.*l2.*m1.*m2.*m3.*t16.*t60.*3.0-g.*l1.*l3.*m1.*m2.*m3.*t14.*t70.*(3.0./4.0)+g.*l1.*l3.*m1.*m2.*m3.*t14.*t71.*(9.0./4.0)+g.*l1.*l3.*m1.*m2.*m3.*t14.*t72.*(3.0./4.0)+dth1.*dth2.*l1.*m1.*t14.*t16.*t19.*t26.*1.2e+1-dth1.*dth2.*l1.*m2.*t14.*t16.*t19.*t25.*1.5e+1-dth1.*dth2.*l1.*m3.*t14.*t16.*t18.*t25.*6.0+dth1.*dth3.*l1.*m1.*t14.*t16.*t19.*t26.*6.0+dth1.*dth3.*l1.*m2.*t14.*t16.*t19.*t26.*9.0+dth2.*dth3.*l1.*m1.*t14.*t16.*t19.*t26.*6.0+dth2.*dth3.*l1.*m2.*t14.*t16.*t19.*t26.*9.0+dth1.*dth2.*l1.*m2.*t14.*t16.*t26.*t52-dth1.*dth2.*l1.*m2.*t14.*t16.*t19.*t68.*3.0-dth1.*dth3.*l1.*m2.*t14.*t16.*t19.*t68.*3.0-dth2.*dth3.*l1.*m2.*t14.*t16.*t19.*t68.*3.0+l1.*l3.*m1.*m2.*m3.*t6.*t9.*t15.*4.0+l1.*l3.*m1.*m2.*m3.*t6.*t10.*t15.*4.0-l2.*m1.*m2.*m3.*t5.*t9.*t12.*t16.*4.0+l3.*m1.*m2.*m3.*t9.*t12.*t14.*t29-l3.*m1.*m2.*m3.*t9.*t12.*t14.*t62.*3.0+dth1.*dth2.*l1.*l3.*m1.*m2.*m3.*t6.*t15.*8.0).*-6.0)./(l1.*m3.*t16);

    end
end