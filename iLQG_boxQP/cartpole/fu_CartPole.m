function fu = fu_CartPole(sys, x, u)

mp = sys.mp;
mc = sys.mc;
l = sys.l;

fu = zeros(size(x, 1), size(u, 1), size(x, 2));
fu(2,1,:) = 1.0./(mc + mp*sin(x(3,:)).^2);
fu(2,2,:) = -cos(x(3,:))./(l*(mc + mp*sin(x(3,:)).^2));
fu(4,1,:) = -cos(x(3,:))./(l*(mc + mp*sin(x(3,:)).^2));
fu(4,2,:) = (1/l*(mc/mp+1))./(l*(mc + mp*sin(x(3,:)).^2));
end
