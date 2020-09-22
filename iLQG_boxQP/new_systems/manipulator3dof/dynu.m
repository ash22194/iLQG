function out = dynu(sys, x, u)
%DYNU
%    OUT = DYNU(SYS, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    21-Sep-2020 10:38:56

m1 = sys.m(1);
m2 = sys.m(2);
m3 = sys.m(3);
l1 = sys.l(1);
l2 = sys.l(2);
l3 = sys.l(3);

th2 = x(2,:);
th3 = x(3,:);

t2 = cos(th2);
t3 = cos(th3);
t4 = sin(th2);
t5 = sin(th3);
t6 = th2+th3;
t7 = l1.^2;
t8 = l2.^2;
t9 = l3.^2;
t10 = m2.^2;
t11 = m3.^2;
t12 = th2.*2.0;
t13 = th3.*2.0;
t18 = 1.0./l1;
t20 = 1.0./l2;
t22 = 1.0./l3;
t23 = -th3;
t24 = l2.*m2.*8.0;
t25 = l2.*m3.*1.5e+1;
t26 = m1.*m2.*1.6e+1;
t27 = m1.*m3.*3.0e+1;
t34 = m2.*m3.*7.5e+1;
t35 = l1.*l3.*m1.*4.8e+1;
t41 = l1.*l3.*m3.*9.0e+1;
t42 = l1.*l3.*m2.*1.44e+2;
t14 = cos(t12);
t15 = cos(t13);
t16 = t3.^2;
t17 = cos(t6);
t19 = 1.0./t7;
t21 = 1.0./t8;
t28 = t6+th3;
t29 = t6+th2;
t32 = t11.*1.8e+1;
t33 = t10.*3.0e+1;
t36 = t23+th2;
t37 = t6.*2.0;
t39 = l1.*m2.*t2.*1.2e+1;
t40 = l1.*m3.*t2.*1.5e+1;
t43 = l3.*m2.*t2.*7.2e+1;
t44 = l3.*m3.*t2.*1.44e+2;
t50 = l1.*l2.*m1.*t3.*7.2e+1;
t52 = l2.*l3.*m3.*t2.*9.0e+1;
t53 = l1.*l2.*m3.*t3.*1.08e+2;
t54 = l1.*l2.*m2.*t3.*1.62e+2;
t55 = l2.*m2.*t2.*t3.*3.6e+1;
t56 = l2.*m2.*t4.*t5.*7.2e+1;
t57 = l2.*m3.*t4.*t5.*2.16e+2;
t74 = l3.*m3.*t3.*t4.*t5.*1.08e+2;
t30 = cos(t28);
t31 = cos(t29);
t38 = cos(t36);
t45 = cos(t37);
t46 = l2.*m3.*t15.*9.0;
t47 = m1.*m3.*t15.*1.8e+1;
t48 = m2.*m3.*t15.*2.7e+1;
t49 = m2.*m3.*t14.*4.5e+1;
t51 = l2.*t43;
t58 = t10.*t14.*1.8e+1;
t59 = t14.*t32;
t65 = m2.*t8.*t17.*1.8e+1;
t67 = t11.*t14.*-1.8e+1;
t70 = m3.*t8.*t17.*1.08e+2;
t71 = l3.*m3.*t2.*t16.*1.08e+2;
t60 = -t46;
t61 = -t47;
t62 = -t48;
t63 = -t49;
t64 = l1.*m3.*t30.*9.0;
t66 = -t58;
t69 = -t65;
t72 = l1.*l2.*m2.*t31.*5.4e+1;
t73 = l2.*l3.*m3.*t30.*5.4e+1;
t75 = -t70;
t76 = -t71;
t77 = m2.*m3.*t45.*9.0;
t80 = l1.*l2.*m3.*t31.*1.08e+2;
t81 = l1.*l3.*m3.*t45.*5.4e+1;
t82 = m2.*t8.*t38.*5.4e+1;
t85 = m3.*t8.*t38.*1.08e+2;
t68 = -t64;
t78 = -t72;
t79 = -t73;
t83 = -t80;
t84 = -t81;
t87 = t43+t44+t55+t56+t57+t74+t76;
t88 = t26+t27+t32+t33+t34+t61+t62+t63+t66+t67+t77;
t86 = t24+t25+t39+t40+t60+t68;
t89 = 1.0./t88;
t92 = t35+t41+t42+t50+t51+t52+t53+t54+t69+t75+t78+t79+t82+t83+t84+t85;
t90 = t19.*t20.*t86.*t89.*6.0;
t93 = t18.*t20.*t22.*t87.*t89;
t94 = t18.*t21.*t22.*t89.*t92;
t91 = -t90;
t95 = -t94;

out = zeros(size(x,1), size(u,1), size(x,2));
out(4,1,:) = t19.*t89.*(m2.*8.0+m3.*1.5e+1-m3.*t15.*9.0).*6.0;
out(5,1,:) = t91;
out(6,1,:) = t93;

out(4,2,:) = t91;
out(5,2,:) = t19.*t21.*t89.*(m1.*t7.*8.0+m2.*t7.*2.4e+1+m2.*t8.*8.0+m3.*t7.*1.5e+1+m3.*t8.*1.5e+1-m3.*t8.*t15.*9.0-m3.*t7.*t45.*9.0+l1.*l2.*m2.*t2.*2.4e+1+l1.*l2.*m3.*t2.*3.0e+1-l1.*l2.*m3.*t30.*1.8e+1).*6.0;
out(6,2,:) = t95;

out(4,3,:) = t93;
out(5,3,:) = t95;
out(6,3,:) = (t21.*t89.*(t8.*t10.*1.5e+1+t8.*t11.*3.6e+1+t9.*t11.*1.5e+1+m1.*m2.*t8.*8.0+m1.*m3.*t8.*2.4e+1+m1.*m3.*t9.*8.0+m2.*m3.*t8.*6.0e+1+m2.*m3.*t9.*2.4e+1-t8.*t10.*t14.*9.0-t8.*t11.*t14.*3.6e+1-t9.*t11.*t45.*9.0+l2.*l3.*t3.*t11.*3.6e+1-l2.*l3.*t11.*t31.*3.6e+1-m2.*m3.*t8.*t14.*3.6e+1+l2.*l3.*m1.*m3.*t3.*2.4e+1+l2.*l3.*m2.*m3.*t3.*5.4e+1-l2.*l3.*m2.*m3.*t31.*1.8e+1).*6.0)./(m3.*t9);

end