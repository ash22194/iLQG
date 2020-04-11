function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = block_dyn_cst(sys, x, u, full_DDP)
    final = isnan(u(1,:));
%     u(:,final) = 0;
    u(:,final) = -9.81/11;
    if (nargout == 2)
        
        f = (eye(size(x,1)) + [0, sys.dt;0, 0])*x + [0; sys.dt/sys.m]*(u + 9.81/11);
        c = diag(x' * sys.Q * x + (u + 9.81/11)' * sys.R * (u + 9.81/11))' * sys.dt;
    else
        
        f = [];
        X_DIM = size(x, 1);
        fx = repmat(eye(X_DIM) + [0, sys.dt;0, 0], [1, 1, size(x,2)]);
        fu = repmat([0; sys.dt/sys.m], [1, 1, size(x,2)]);

        if (full_DDP)
            [fxx,fxu,fuu] = deal([]);
        else
            [fxx,fxu,fuu] = deal([]);
        end

        c = [];
        cx = 2 * sys.Q * x * sys.dt;
        cu = 2 * sys.R * (u + 9.81/11) * sys.dt;
        cxx = repmat(2 * sys.Q, [1, 1, size(x,2)]) * sys.dt;
        cuu = repmat(2 * sys.R, [1, 1, size(x,2)]) * sys.dt;
        cxu = permute(zeros(size(u,1), size(x,1), size(x,2)), [2, 1, 3]) * sys.dt;
    end
end