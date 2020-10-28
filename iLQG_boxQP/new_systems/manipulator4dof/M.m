function out = M(sys, x, u)
%M
%    OUT = M(SYS, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    19-Oct-2020 18:27:49

l1 = sys.l(1);
l2 = sys.l(2);
l3 = sys.l(3);
l4 = sys.l(4);

m1 = sys.m(1);
m2 = sys.m(2);
m3 = sys.m(3);
m4 = sys.m(4);

th2 = x(2,:);
th3 = x(3,:);
th4 = x(4,:);

t2 = cos(th2);
t3 = cos(th3);
t4 = cos(th4);
t5 = th2+th3;
t6 = th3+th4;
t7 = l4.*2.0;
t8 = l1.^2;
t9 = l2.^2;
t10 = l3.^2;
t11 = l4.^2;
t12 = cos(t5);
t13 = cos(t6);
t14 = t5+th4;
t16 = m3.*t9;
t17 = m4.*t9;
t18 = m4.*t10;
t19 = l3.*t4.*3.0;
t20 = l1.*l2.*m3.*t2;
t21 = l1.*l2.*m4.*t2;
t22 = l2.*l3.*m3.*t3;
t23 = l2.*l3.*m4.*t3;
t24 = l3.*l4.*m4.*t4;
t29 = (m2.*t9)./3.0;
t30 = (m3.*t10)./3.0;
t31 = (m4.*t11)./3.0;
t33 = (l1.*l2.*m2.*t2)./2.0;
t15 = cos(t14);
t25 = l2.*t13.*3.0;
t26 = t23.*2.0;
t27 = l1.*l3.*m4.*t12;
t28 = l2.*l4.*m4.*t13;
t34 = t22./2.0;
t35 = t7+t19;
t36 = (l1.*l3.*m3.*t12)./2.0;
t32 = l1.*t15.*3.0;
t37 = t28./2.0;
t38 = (l1.*l4.*m4.*t15)./2.0;
t39 = (l4.*m4.*t35)./6.0;
t40 = t25+t35;
t41 = (l4.*m4.*t40)./6.0;
t42 = t32+t40;
t44 = t18+t23+t24+t30+t31+t34+t37;
t46 = t16+t17+t18+t20+t21+t22+t24+t26+t27+t28+t29+t30+t31+t33+t36+t38;
t43 = (l4.*m4.*t42)./6.0;
t45 = t27+t36+t38+t44;

out = zeros(4,4,size(x,2));
out(1,1,:) = t16+t17+t18+t20.*2.0+t21.*2.0+t22+t24+t26+t27.*2.0+t28+t29+t30+t31+(m1.*t8)./3.0+m2.*t8+m3.*t8+m4.*t8+l1.*l2.*m2.*t2+l1.*l3.*m3.*t12+l1.*l4.*m4.*t15;
out(2,1,:) = t46;
out(3,1,:) = t45;
out(4,1,:) = t43;
out(1,2,:) = t46;
out(1,3,:) = t45;
out(1,4,:) = t43;

out(2,2,:) = t16+t17+t18+t22+t24+t26+t28+t29+t30+t31;
out(3,2,:) = t44;
out(4,2,:) = t41;
out(2,3,:) = t44;
out(2,4,:) = t41;

out(3,3,:) = t18+t24+t30+t31;
out(4,3,:) = t39;
out(3,4,:) = t39;

out(4,4,:) = t31;
end
