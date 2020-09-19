function out = dyn(sys, x, u)
%DYN
%    OUT = DYN(SYS,X,U)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Sep-2020 19:09:32

mc = sys.mc;
mp = sys.mp;
l = sys.l;
g = sys.g;

x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);

u1 = u(1,:);
u2 = u(2,:);

t2 = cos(x3);
t3 = sin(x3);
t4 = x4.^2;
t6 = 1.0./l;
t5 = t3.^2;
t7 = mp.*t5;
t8 = mc+t7;
t9 = 1.0./t8;
out = [x2;
       t6.*t9.*(l.*u1-t2.*u2+mp.*t3.*t4.*1.0./t6.^2+g.*l.*mp.*t2.*t3);
       x4;
       -t6.*t9.*(t2.*u1-t6.*u2.*(mc./mp+1.0))-t6.*t9.*(g.*t3.*(mc+mp)+l.*mp.*t2.*t3.*t4)];

end
