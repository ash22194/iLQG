function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = cartpole_dyn_second_wraparound_cst(sys, x, u, k, K, xn, full_DDP)
    
    final = isnan(u(1,:));
    u(:,final) = 0;
    if (nargout == 2)
        f = f_CartPoleSecond_WrapAround_finite(sys, x, u, k, K, xn, sys.dt);
        c = l_CartPoleSecond_WrapAround(sys, x, u, k, K, xn) * sys.dt;
    else
        
        f = [];
        X_DIM = size(x, 1);
        fx = repmat(eye(X_DIM), [1,1,size(x,2)]) ...
              + fx_CartPoleSecond_WrapAround(sys, x, u, k, K, xn) * sys.dt;
        fu = fu_CartPoleSecond(sys, x, u, k, K, xn) * sys.dt;
        if (full_DDP)
            [fxx,fxu,fuu] = deal([]);
        else
            [fxx,fxu,fuu] = deal([]);
        end

        c = [];
        cx = dldx_CartPoleSecond_WrapAround(sys, x, u, k, K, xn) * sys.dt;
        cu = dldu_CartPoleSecond_WrapAround(sys, x, u, k, K, xn) * sys.dt;
        cxx = ddldx_CartPoleSecond(sys, x, u, k, K, xn) * sys.dt;
        cuu = ddldu_CartPoleSecond(sys, x, u, k, K, xn) * sys.dt;
        cxu = permute(ddldudx_CartPoleSecond(sys, x, u, k, K, xn), [2, 1, 3]) * sys.dt;
    end
end