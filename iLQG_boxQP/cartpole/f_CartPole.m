function dx = f_CartPole(sys, x, u)

mp = sys.mp;
mc = sys.mc;
l = sys.l;
g = sys.g;

zero_dyn = [x(2,:);
            (mp*l*(x(4,:).^2).*sin(x(3,:)) + 0.5*mp*g*sin(2*x(3,:)))./(mc + mp*sin(x(3,:)).^2);
            x(4,:);
            ((-0.5*mp*l*x(4,:).^2).*sin(2*x(3,:)) - (mc+mp)*g*sin(x(3,:)))./(l*(mc + mp*sin(x(3,:)).^2))];
act_dyn = [zeros(1, size(x,2));
           (u(1,:) - u(2,:).*cos(x(3,:))/l)./(mc + mp*sin(x(3,:)).^2);
           zeros(1, size(x,2));
           (u(2,:)/l*(mc/mp+1) - u(1,:).*cos(x(3,:)))./(l*(mc + mp*sin(x(3,:)).^2))];

dx = zero_dyn + act_dyn;

end
