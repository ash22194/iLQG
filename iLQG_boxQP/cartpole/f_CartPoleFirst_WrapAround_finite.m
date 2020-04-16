function x_new = f_CartPoleFirst_WrapAround_finite(sys, x, u, dt)

wraparound_dim = find(sys.X_DIMS_FREE==3);
x(wraparound_dim, :) = mod(x(wraparound_dim,:), 2*pi);

k1 = f_CartPoleFirst(sys, x, u);
k2 = f_CartPoleFirst(sys, x + 0.5 * dt * k1, u);
k3 = f_CartPoleFirst(sys, x + 0.5 * dt * k2, u);
k4 = f_CartPoleFirst(sys, x + dt * k3, u);

x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

end

