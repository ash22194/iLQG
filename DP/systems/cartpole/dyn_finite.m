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
        grid_size = ones(length(x),1);
        if (length(x) == (sys.X_DIMS+1))
            grid_size([sys.X_DIMS_FREE; (sys.X_DIMS+1)]) = size(x{sys.X_DIMS_FREE(1)});
        else
            grid_size(sys.X_DIMS_FREE) = size(x{sys.X_DIMS_FREE(1)});
        end
        
        grid_size = int32(grid_size);
    else
        grid_size = int32(sys.grid_size);
    end
    
    limits = sys.limits;
    
    out = cell(4,1);
    [out{:}] = dyn_mex_finite(x{:}, u{:}, sys.mc, sys.mp, ...
                              sys.l, sys.g, dt, ...
                              grid_size, active_actions, limits, sys.X_DIMS_FREE(:));

    if (length(x)==(sys.X_DIMS+1))
        out = cat(1, out, {min(x{sys.X_DIMS+1} + dt, sys.limits_t(2))});
    end
end