function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = cartpole_dyn_cst(sys, x, u, full_DDP)

    if (nargout == 2)
        final = isnan(u(1,:));
        u(:,final) = 0;
        f = f_CartPole_finite(sys, x, u, sys.dt);
        c = l_CartPole(sys, x, u) * sys.dt;
%         f = f_CartPole_WrapAround_finite(sys, x, u, sys.dt);
%         c = l_CartPole_WrapAround(sys, x, u) * sys.dt;
    else
        
        f = [];
        X_DIM = size(x, 1);
        fx = repmat(eye(X_DIM), [1,1,size(x,2)]) + fx_CartPole(sys, x, u) * sys.dt;
        fu = fu_CartPole(sys, x, u) * sys.dt;
%         xu_dyn  = @(xu) f_CartPole_WrapAround_finite(sys, xu(1:4,:), xu(5:6,:), sys.dt);
%         J       = finite_difference(xu_dyn, [x; u]);
%         fx      = J(:,1:4,:);
%         fu      = J(:,5:6,:);
        if (full_DDP)
            [fxx,fxu,fuu] = deal([]);
        else
            [fxx,fxu,fuu] = deal([]);
        end

        c = [];
        cx = dldx_CartPole(sys, x, u) * sys.dt;
%         cx = dldx_CartPole_WrapAround(sys, x, u) * sys.dt;
        cu = dldu_CartPole(sys, x, u) * sys.dt;
        cxx = ddldx_CartPole(sys, x, u) * sys.dt;
        cuu = ddldu_CartPole(sys, x, u) * sys.dt;
        cxu = permute(ddldudx_CartPole(sys, x, u), [2, 1, 3]) * sys.dt;
    end
end

function J = finite_difference(fun, x, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 3
    h = 2^-17;
end

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));
Y       = fun(X);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);
end

function c = pp(a,b)
c = bsxfun(@plus,a,b);
end