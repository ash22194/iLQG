function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = cartpole_dyn_second_wraparound_cst(sys, x, u, k, K, xn, full_DDP)
    
    final = isnan(u(1,:));
    u(:,final) = 0;
    if (nargout == 2)
        f = f_CartPoleSecond_WrapAround_finite(sys, x, u, k, K, xn, sys.dt);
        c = l_CartPoleSecond_WrapAround(sys, x, u, k, K, xn) * sys.dt;
    else
        
        f = [];
        X_DIMS_FREE = sys.X_DIMS_FREE;
        xu_dyn  = @(xu, k, K, xn) f_CartPoleSecond_WrapAround_finite(sys, xu(1:length(X_DIMS_FREE),:), ...
                                                                xu((1+length(X_DIMS_FREE)):end,:), k, K, xn, sys.dt);
        J       = finite_difference(xu_dyn, [x; u], k, K, xn);
        fx      = J(:,1:length(X_DIMS_FREE),:);
        fu      = J(:,(1+length(X_DIMS_FREE)):end,:);
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

function J = finite_difference(fun, x, kF, KF, xnF, h)
% simple finite-difference derivatives
% assumes the function fun() is vectorized

if nargin < 6
    h = 2^-17;
end

[n, K]  = size(x);
H       = [zeros(n,1) h*eye(n)];
H       = permute(H, [1 3 2]);
X       = pp(x, H);
X       = reshape(X, n, K*(n+1));

dimkF   =  size(kF,1);
HkF     = zeros(dimkF, n+1);
HkF     = permute(HkF, [1 3 2]);
kF      = pp(kF, HkF);
kF      = reshape(kF, dimkF, K*(n+1));

dimxnF  = size(xnF,1);
HxnF    = zeros(dimxnF, n+1);
HxnF    = permute(HxnF, [1 3 2]);
xnF     = pp(xnF, HxnF);
xnF     = reshape(xnF, dimxnF, K*(n+1));

dimKF1  =  size(KF,1);
dimKF2  =  size(KF,2);
HKF     = zeros(dimKF1, dimKF2, n+1);
HKF     = permute(HKF, [1 2 4 3]);
KF      = pp(KF, HKF);
KF      = reshape(KF, dimKF1, dimKF2, K*(n+1));

Y       = fun(X, kF, KF, xnF);
m       = numel(Y)/(K*(n+1));
Y       = reshape(Y, m, K, n+1);
J       = pp(Y(:,:,2:end), -Y(:,:,1)) / h;
J       = permute(J, [1 3 2]);
end

function c = pp(a,b)
c = bsxfun(@plus,a,b);
end