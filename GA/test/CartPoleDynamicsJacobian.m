function out1 = CartPoleDynamicsJacobian(x, u)
%CARTPOLEDYNAMICSJACOBIAN
%    OUT1 = CARTPOLEDYNAMICSJACOBIAN(U1,U2,X3,X4)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    02-Sep-2020 11:55:35

x3 = x(3);
x4 = x(4);

u1 = u(1);
u2 = u(2);

t2 = cos(x3);
t3 = sin(x3);
t4 = x3.*2.0;
t5 = x4.^2;
t6 = cos(t4);
t7 = sin(t4);
t8 = t3.^2;
t9 = t8+5.0;
t12 = t8.*(9.0./1.0e+1);
t10 = 1.0./t9;
t13 = t12+9.0./2.0;
t11 = t10.^2;
t14 = 1.0./t13;
t15 = t14.^2;
out1 = reshape([0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,t10.*(t6.*(9.81e+2./1.0e+2)+t2.*t5.*(9.0./1.0e+1))+t3.*t10.*u2.*(1.0e+1./9.0)-t2.*t3.*t11.*(t7.*(9.81e+2./2.0e+2)+t3.*t5.*(9.0./1.0e+1)).*2.0-t2.*t3.*t11.*(u1-t2.*u2.*(1.0e+1./9.0)).*2.0,0.0,-t14.*(t2.*5.886e+1+t5.*t6.*(9.0./1.0e+1))+t3.*t14.*u1+t2.*t3.*t15.*(t3.*5.886e+1+t5.*t7.*(9.0./2.0e+1)).*(9.0./5.0)-t2.*t3.*t15.*(u2.*(2.0e+1./3.0)-t2.*u1).*(9.0./5.0),0.0,t3.*t10.*x4.*(9.0./5.0),1.0,t7.*t14.*x4.*(-9.0./1.0e+1),0.0,t10,0.0,-t2.*t14,0.0,t2.*t10.*(-1.0e+1./9.0),0.0,t14.*(2.0e+1./3.0)],[4,6]);
end
