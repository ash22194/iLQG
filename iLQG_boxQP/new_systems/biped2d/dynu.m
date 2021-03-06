function out = dynu(sys, x, u)
%DYNU
%    OUT = DYNU(SYS, X, U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    15-Sep-2020 22:32:10

m = sys.m;
I = sys.I;
d = sys.d;
df = sys.df;
l0 = sys.l0;

x1 = x(1, :);
x2 = x(2, :);
x5 = x(5, :);

t2 = cos(x2);
t3 = sin(x2);
l2 = sqrt((x1.*t2 + df).^2 + (x1.*t3).^2);
contact1 = x1 <= l0;
contact2 = l2 <= l0;

t4 = df.^2;
t5 = x1.^2;
t7 = 1.0./I;
t9 = 1.0./m;
t10 = -x5;
t11 = 1.0./x1;
t6 = t2.^2;
t8 = t2.*x1;
t14 = t10+x2;
t12 = df.*t8.*2.0;
t13 = df+t8;
t15 = t6-1.0;
t16 = t4+t5+t12;
t17 = 1.0./t16;
t18 = 1.0./sqrt(t16);
t19 = t13.*t18;
t21 = t5.*t15.*t17;
t20 = acos(t19);
t23 = -t21;
t22 = -t20;
t25 = sqrt(t23);
t24 = t22+x5;

t2 = t2.*contact1;
t3 = t3.*contact1;

t17 = t17.*contact2;
t19 = t19.*contact2;
t25 = t25.*contact2;

out = zeros(size(x, 1), size(u, 1), size(u, 2));
out(3,1,:) = t2.*t9;
out(4,1,:) = t3.*t9;
out(6,1,:) = d.*t7.*cos(t14).*contact1;

out(3,2,:) = t9.*t19;
out(4,2,:) = t9.*t25;
out(6,2,:) = d.*t7.*cos(t24).*contact2;

out(3,3,:) = t3.*t9.*t11;
out(4,3,:) = -t2.*t9.*t11;
out(6,3,:) = t7.*t11.*(x1+d.*sin(t14)).*contact1;

out(3,4,:) = t9.*t18.*t25;
out(4,4,:) = -t9.*t13.*t17;
out(6,4,:) = -t7.*t18.*(d.*sin(t24)-1.0./t18).*contact2;

end
