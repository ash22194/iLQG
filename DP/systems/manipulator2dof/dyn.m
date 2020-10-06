function out = dyn(sys, x, u)
%DYN
%    OUT = DYN(SYS, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    21-Sep-2020 10:14:09

m1 = sys.m(1);
m2 = sys.m(2);
l1 = sys.l(1);
l2 = sys.l(2);
g = sys.g;
X_DIMS_FREE = sys.X_DIMS_FREE; % Assuming they are sorted

th1 = x{1};
th2 = x{2};
dth1 = x{3};
dth2 = x{4};

u1 = u{1};
u2 = u{2};

t2 = cos(th1);
t3 = cos(th2);
t4 = sin(th1);
t5 = sin(th2);
t6 = dth1.^2;
t7 = dth2.^2;
t8 = l1.^2;
t9 = l1.^3;
t10 = l2.^2;
t11 = l2.^3;
t12 = m2.^2;
t13 = th2.*2.0;
t14 = cos(t13);
t15 = sin(t13);
t16 = 1.0./t8;

out = cell(length(X_DIMS_FREE));

for xxi=1:1:length(X_DIMS_FREE)
    xx = X_DIMS_FREE(xxi);
    if (xx==1)
        out{xxi} = dth1;
    elseif (xx==2)
        out{xxi} = dth2;
    elseif (xx==3)
        out{xxi} = (t16.*(l2.*u1.*4.0-l2.*u2.*4.0-l1.*t3.*u2.*6.0-g.*l1.*l2.*m1.*t4.*2.0-g.*l1.*l2.*m2.*t4.*(5.0./2.0)+l1.*m2.*t5.*t6.*t10.*2.0+l1.*m2.*t5.*t7.*t10.*2.0+l2.*m2.*t6.*t8.*t15.*(3.0./2.0)+g.*l1.*l2.*m2.*sin(t13+th1).*(3.0./2.0)+dth1.*dth2.*l1.*m2.*t5.*t10.*4.0).*3.0)./(l2.*(m1.*4.0+m2.*1.2e+1-m2.*t3.^2.*9.0));
    elseif (xx==4)
        out{xxi} = (t16.*(m1.*t8.*u2.*-8.0-m2.*t8.*u2.*2.4e+1+m2.*t10.*u1.*8.0-m2.*t10.*u2.*8.0+l1.*l2.*m2.*t3.*u1.*1.2e+1-l1.*l2.*m2.*t3.*u2.*2.4e+1-g.*l1.*t4.*t10.*t12.*5.0+l2.*t5.*t6.*t9.*t12.*1.2e+1+l1.*t5.*t6.*t11.*t12.*4.0+l1.*t5.*t7.*t11.*t12.*4.0+t6.*t8.*t10.*t12.*t15.*6.0+t7.*t8.*t10.*t12.*t15.*3.0+l2.*m1.*m2.*t5.*t6.*t9.*4.0+g.*l2.*t2.*t5.*t8.*t12.*1.2e+1+g.*l1.*t2.*t10.*t12.*t15.*3.0+g.*l1.*t4.*t10.*t12.*t14.*3.0+dth1.*dth2.*l1.*t5.*t11.*t12.*8.0-g.*l1.*m1.*m2.*t4.*t10.*4.0+dth1.*dth2.*t8.*t10.*t12.*t15.*6.0+g.*l2.*m1.*m2.*t2.*t5.*t8.*4.0-g.*l2.*m1.*m2.*t3.*t4.*t8.*2.0).*-3.0)./(m2.*t10.*(m1.*8.0+m2.*1.5e+1-m2.*t14.*9.0));
    end
end

end