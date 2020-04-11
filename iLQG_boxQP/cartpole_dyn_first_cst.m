function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = cartpole_dyn_first_cst(sys, x, u, full_DDP)
    
    final = isnan(u(1,:));
    u(:,final) = 0;
    if (nargout == 2)
        f = f_CartPoleFirst_finite(sys, x, u, sys.dt);
        c = l_CartPoleFirst(sys, x, u) * sys.dt;
    else
        
        f = [];
        X_DIM = size(x, 1);
        fx = repmat(eye(X_DIM), [1,1,size(x,2)]) + fx_CartPoleFirst(sys, x, u) * sys.dt;
        fu = fu_CartPoleFirst(sys, x, u) * sys.dt;
        if (full_DDP)
            [fxx,fxu,fuu] = deal([]);
        else
            [fxx,fxu,fuu] = deal([]);
        end

        c = [];
        cx = dldx_CartPoleFirst(sys, x, u) * sys.dt;
        cu = dldu_CartPoleFirst(sys, x, u) * sys.dt;
        cxx = ddldx_CartPoleFirst(sys, x, u) * sys.dt;
        cuu = ddldu_CartPoleFirst(sys, x, u) * sys.dt;
        cxu = permute(ddldudx_CartPoleFirst(sys, x, u), [2, 1, 3]) * sys.dt;
    end
end