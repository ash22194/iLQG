function out = dyn_finite_rk4_gpumex(sys, x, u, dt)

X_DIMS_FREE = sys.X_DIMS_FREE;
limits      = sys.limits;
M = num2cell(sys.m);
L = num2cell(sys.l);
g = sys.g;
grid_size = int32(size(x{1}));

k1 = cell(length(x), 1);
[k1{:}] = dyn_mex(x{:}, u{:}, M{:}, L{:}, g, grid_size);
q = x;
for xxi=1:1:length(X_DIMS_FREE)
    xx = X_DIMS_FREE(xxi);
%     q{xx} = q{xx} + 0.5 * dt * k1{xxi};
    q{xx} = arrayfun(@step, q{xx}, k1{xxi}, 0.5*dt);
end

k2 = cell(length(x),1);
[k2{:}] = dyn_mex(q{:}, u{:}, M{:}, L{:}, g, grid_size);
q = x;
for xxi=1:1:length(X_DIMS_FREE)
    xx = X_DIMS_FREE(xxi);
%     q{xx} = q{xx} + 0.5 * dt * k2{xxi};
    q{xx} = arrayfun(@step, q{xx}, k2{xxi}, 0.5*dt);
end

k3 = cell(length(x),1);
[k3{:}] = dyn_mex(q{:}, u{:}, M{:}, L{:}, g, grid_size);
q = x;
for xxi=1:1:length(X_DIMS_FREE)
    xx = X_DIMS_FREE(xxi);
%     q{xx} = q{xx} + dt * k3{xxi};
    q{xx} = arrayfun(@step, q{xx}, k3{xxi}, dt);
end

k4 = cell(length(x),1);
[k4{:}] = dyn_mex(q{:}, u{:}, M{:}, L{:}, g, grid_size);
out = x;

for xxi=1:1:length(X_DIMS_FREE)
    xx = X_DIMS_FREE(xxi);
%     out{xx} = out{xx} + dt / 6.0 * (k1{xxi} + 2.0 * k2{xxi} + 2.0 * k3{xxi} + k4{xxi});
    out{xx} = arrayfun(@rk4, out{xx}, k1{xxi}, k2{xxi}, k3{xxi}, k4{xxi}, dt);
    
    % Enforce state bounds!!
%     out{xx} = max(limits(xx, 1), min(limits(xx, 2), out{xx}));
    out{xx} = arrayfun(@bound, out{xx}, limits(xx, 1), limits(xx, 2));
end

    function y_ = step(y, k, delta)
        y_ = y + delta * k;
    end
    
    function y_ = rk4(y, x1, x2, x3, x4, delta)
        y_ = y + delta / 6.0 * (x1 + 2.0 * x2 + 2.0 * x3 + x4);
    end
    
    function y_ = bound(y, mi, ma)
        y_ = max(mi, min(ma, y));
    end
end
