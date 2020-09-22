clear;
close all;
clc;

%% 

n = 2;
addpath(sprintf('../iLQG_boxQP/new_systems/manipulator%ddof', n));
load('sys.mat');

assert(isfield(sys, 'name') && strcmp(sys.name, sprintf('manipulator%ddof', n)), 'Check sys.mat');

sys.m = [5; 1]; % kg
sys.l = [0.24; 0.9]; % m
sys.g = 9.81; % m/s^2
sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(1) = pi;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.dt = 0.001;

sys.fxu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
fxu = sys.fxu_func(sys.l_point, sys.u0);
sys.A = fxu(:,1:sys.X_DIMS); 
sys.B = fxu(:,(sys.X_DIMS+1):(sys.X_DIMS + sys.U_DIMS));

sys.gamma_ = 0.997;
sys.lambda_ = (1 - sys.gamma_) / sys.dt;
sys.Q = diag([25*ones(1, sys.n), 0.02*ones(1, sys.n)]);
sys.R = diag(0.001*ones(1, sys.n));

[K_joint, S_joint, ~] = lqr(sys.A - eye(size(sys.A,1))*sys.lambda_/2, sys.B, sys.Q, sys.R, zeros(size(sys.A,1), size(sys.B,2)));
sys.S =  sym('S', [sys.X_DIMS, sys.X_DIMS]);
sys.S = tril(sys.S,0) + tril(sys.S,-1).';
sys.a = sym('a', [sys.X_DIMS, 1]);
sys.b = sym('b', [sys.X_DIMS, 1]);
sys.err_lqr = sys.x.'*sys.S*sys.x - sys.x.'*S_joint*sys.x;
for ii=1:1:sys.X_DIMS
    sys.err_lqr = int(sys.err_lqr, sys.x(ii), [sys.a(ii), sys.b(ii)]);
end
% matlabFunction(simplify(sys.err_lqr, 'IgnoreAnalyticConstraints', true), 'File', sprintf('Manipulator%dDOFErrLQR.m', n));
sys.err_lqr_func = @(S,a,b) Manipulator2DOFErrLQR(S,a,b);

sys.numPoints = 31*ones(1, sys.X_DIMS);
sys.limits = [zeros(n, 1), 2*pi*ones(n, 1);
              -3*ones(n, 1), 3*ones(n, 1)];
sys.state_bounds = [-pi/3*ones(n, 1) + sys.l_point(1:n), pi/3*ones(n, 1) + sys.l_point(1:n);
                    -ones(n, 1), ones(n, 1)];
sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1));

p = [0, 1;
     0, 2];
s = [1, 0, 1, 0;
     0, 1, 0, 1];

err_lqr = computeLQRMeasure(sys, p, s);
err_compute = computeComplexityEstimates(sys, p, s);
