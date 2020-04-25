function x_new = f_Biped2DSecond_finite(sys, x, u, k, K, xn, dt)

k1 = f_Biped2DSecond(sys, x, u, k, K, xn);
k2 = f_Biped2DSecond(sys, x + 0.5 * dt * k1, u, k, K, xn);
k3 = f_Biped2DSecond(sys, x + 0.5 * dt * k2, u, k, K, xn);
k4 = f_Biped2DSecond(sys, x + dt * k3, u, k, K, xn);

x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

end
