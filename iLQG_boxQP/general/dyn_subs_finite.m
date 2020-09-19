function x_new = dyn_subs_finite(sys, x, u, sub_policies, dt)
    
k1 = dyn_subs(sys, x, u, sub_policies);
k2 = dyn_subs(sys, x + 0.5 * dt * k1, u, sub_policies);
k3 = dyn_subs(sys, x + 0.5 * dt * k2, u, sub_policies);
k4 = dyn_subs(sys, x + dt * k3, u, sub_policies);

x_new = x + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

end