function out = dyn_finite(sys, x, u, dt)

    if (~isfield(sys, 'active_actions'))
        active_actions = zeros(sys.U_DIMS,1);
        active_actions(sys.U_DIMS_FREE) = 1;
        active_actions(sys.U_DIMS_CONTROLLED) = 1;
        active_actions = int32(active_actions);
    else
        active_actions = int32(sys.active_actions);
    end

    if (~isfield(sys, 'grid_size'))
        grid_size = ones(sys.X_DIMS,1);
        grid_size(sys.X_DIMS_FREE) = size(x{sys.X_DIMS_FREE(1)});
        grid_size = int32(grid_size);
    else
        grid_size = int32(sys.grid_size);
    end

    limits = sys.limits;

    out = cell(8,1);
    [out{:}] = dyn_mex_finite(x{:}, u{:}, sys.mw, sys.Iw(1,1), sys.Iw(2,2), sys.Iw(3,3), ...
                              sys.mf, sys.If(1,1), sys.If(2,2), sys.If(3,3), ...
                              sys.mt, sys.It(1,1), sys.It(2,2), sys.It(3,3), ...
                              sys.nw, sys.nt, sys.rw, sys.rf, sys.rt, sys.fcoeff, sys.g, dt, ...
                              grid_size, active_actions, limits, sys.X_DIMS_FREE(:));
end