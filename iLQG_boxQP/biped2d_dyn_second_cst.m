function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = biped2d_dyn_second_cst(sys, x, u, k, K, xn, full_DDP)
    
    u0FREE = sys.u0(sys.U_DIMS_FREE,1);
    u0FIXED = sys.u0(sys.U_DIMS_FIXED,1);
    final = isnan(u(1,:));
    u(:,final) = repmat(u0FREE, [1, sum(final)]);
    if (nargout == 2)
        f = f_Biped2DSecond_finite(sys, x, u, k, K, xn, sys.dt);
        c = l_Biped2DSecond(sys, x, u - u0FREE, k - u0FIXED, K, xn) * sys.dt;
    
    else
        f = [];
        X_DIM = size(x, 1);
        fx = repmat(eye(X_DIM), [1,1,size(x,2)]) ...
              + fx_Biped2DSecond(sys, x, u, k, K, xn) * sys.dt;
        fu = fu_Biped2DSecond(sys, x, u, k, K, xn) * sys.dt;
        if (full_DDP)
            [fxx,fxu,fuu] = deal([]);
        else
            [fxx,fxu,fuu] = deal([]);
        end

        c = [];
        cx = dldx_Biped2DSecond(sys, x, u - u0FREE, k - u0FIXED, K, xn) * sys.dt;
        cu = dldu_Biped2DSecond(sys, x, u - u0FREE, k - u0FIXED, K, xn) * sys.dt;
        cxx = ddldx_Biped2DSecond(sys, x, u - u0FREE, k - u0FIXED, K, xn) * sys.dt;
        cuu = ddldu_Biped2DSecond(sys, x, u - u0FREE, k - u0FIXED, K, xn) * sys.dt;
        cxu = permute(ddldudx_Biped2DSecond(sys, x, u - u0FREE, k - u0FIXED, K, xn), [2, 1, 3]) * sys.dt;
    end
end