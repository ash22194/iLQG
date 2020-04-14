function x_new = f_CartPoleSecond_WrapAround_finite(sys, x, u, k, K, xn, dt)

while (any(x(3,:) > 2*pi))
    gthan = x(3,:) > 2*pi;
    x(3,gthan) = x(3,gthan) - 2*pi;
end
while (any(x(3,:) < 0))
    lthan = x(3,:) < 0;
    x(3,lthan) = x(3,lthan) + 2*pi;
end

k1 = f_CartPoleSecond(sys, x, u, k, K, xn);
k2 = f_CartPoleSecond(sys, x + 0.5 * dt * k1, u, k, K, xn);
k3 = f_CartPoleSecond(sys, x + 0.5 * dt * k2, u, k, K, xn);
k4 = f_CartPoleSecond(sys, x + dt * k3, u, k, K, xn);

x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

end

