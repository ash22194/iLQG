function out = dynxu(sys, x, u)
%DYNXU
%    OUT = DYNXU(SYS, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    26-Sep-2020 16:08:41

m1 = sys.m(1);
m2 = sys.m(2);
l1 = sys.l(1);
l2 = sys.l(2);

th2 = x(2,:);

t2 = cos(th2);
t3 = sin(th2);
t4 = m1.*4.0;
t5 = th2.*2.0;
t8 = 1.0./l1.^2;
t9 = 1.0./l2;
t10 = m2.*1.2e+1;
t6 = t2.^2;
t7 = l1.*t4;
t11 = l1.*t10;
t12 = l2.*t2.*t10;
t13 = m2.*t6.*9.0;
t14 = -t13;
t15 = l1.*t13;
t16 = t4+t10+t14;
t18 = t7+t11+t12+t15;
t17 = 1.0./t16.^2;
t19 = t3.*t8.*t9.*t17.*t18.*1.8e+1;

out = zeros(size(x,1), size(x,1), size(u,1), size(x,2));

out(3,2,1,:) = m2.*t8.*sin(t5).*1.0./(m1.*8.0+m2.*1.5e+1-m2.*cos(t5).*9.0).^2.*-4.32e+2;
out(4,2,1,:) = t19;

out(3,2,2,:) = t19;
out(4,2,2,:) = t3.*t8.*t9.^2.*t17.*(l2.*2.0+l1.*t2.*3.0).*(l1.*m1.*2.0+l1.*m2.*6.0+l2.*m2.*t2.*3.0).*-3.6e+1;

end
