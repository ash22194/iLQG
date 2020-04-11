function x_new = f_CartPoleFirst_finite(sys, x, u, dt)

k1 = f_CartPoleFirst(sys, x, u);
k2 = f_CartPoleFirst(sys, x + 0.5 * dt * k1, u);
k3 = f_CartPoleFirst(sys, x + 0.5 * dt * k2, u);
k4 = f_CartPoleFirst(sys, x + dt * k3, u);

x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

end
