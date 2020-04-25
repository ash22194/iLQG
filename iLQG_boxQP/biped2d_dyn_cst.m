function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = biped2d_dyn_cst(sys, x, u, full_DDP)

    if (nargout == 2)
        final = isnan(u(1,:));
        u(:,final) = 0;
        f = f_Biped2D_finite(sys, x, u, sys.dt);
        c = l_Biped2D(sys, x, u) * sys.dt;

    else

        f = [];
        X_DIM = size(x, 1);
        fx = repmat(eye(X_DIM), [1,1,size(x,2)]) + fx_Biped2D(sys, x, u) * sys.dt;
        fu = fu_Biped2D(sys, x, u) * sys.dt;

        if (full_DDP)
            [fxx,fxu,fuu] = deal([]);
        else
            [fxx,fxu,fuu] = deal([]);
        end

        c = [];
        cx = dldx_Biped2D(sys, x, u) * sys.dt;
        cu = dldu_Biped2D(sys, x, u) * sys.dt;
        cxx = ddldx_Biped2D(sys, x, u) * sys.dt;
        cuu = ddldu_Biped2D(sys, x, u) * sys.dt;
        cxu = permute(ddldudx_Biped2D(sys, x, u), [2, 1, 3]) * sys.dt;
    end
end