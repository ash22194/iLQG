clear;
close all;
clc;

%% 

restoredefaultpath;
system_name = 'cartpole';
addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
load(strcat('../iLQG_boxQP/new_systems/', system_name, '/sys.mat'), 'sys');

sys.X_DIMS = 4; % [x, dx, th, dth]
sys.U_DIMS = 2; % [F, tau]
sys.mc = 5; % kg
sys.mp = 1; % kg
sys.l = 0.9; % m
sys.Q = diag([25, 0.02, 25, 0.02]);
sys.R = diag(0.001*ones(1,2));
sys.g = 9.81; % m/s^2
sys.dt = 0.001;
sys.gamma_ = 0.997;
sys.lambda_ = (1 - sys.gamma_) / sys.dt;
sys.lims = 9*[-ones(2,1), ones(2,1)]; % action limits

% Define decompositions to test
u_x = [];

% Cascaded
p = [0, 1;1, 1];
s = [1, 1, 0, 0;0, 0, 1, 1];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

p = [0, 1;1, 1];
s = [0, 0, 1, 1;1, 1, 0, 0];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

p = [0, 1;1, 1];
s = [0, 0, 0, 0;1, 1, 1, 1];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

p = [2, 1;0, 1];
s = [1, 1, 0, 0;0, 0, 1, 1];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

p = [2, 1;0, 1];
s = [0, 0, 1, 1;1, 1, 0, 0];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

p = [2, 1;0, 1];
s = [1, 1, 1, 1;0, 0, 0, 0];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

% Decoupled
p = [0, 1;0, 2];
s = [1, 1, 0, 0;0, 0, 1, 1];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

p = [0, 1;0, 2];
s = [0, 0, 1, 1;1, 1, 0, 0];
u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];

sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(3) = pi;
sys.goal = sys.l_point;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.fxfu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];

fxfu = sys.fxfu_func(sys.l_point, sys.u0);
sys.A = fxfu(:,1:sys.X_DIMS);
sys.B = fxfu(:,(1+sys.X_DIMS):end);
[~, S_joint, ~] = lqr(sys.A - eye(size(sys.A,1))*sys.lambda_/2, sys.B, sys.Q, sys.R, zeros(size(sys.A,1), size(sys.B,2)));
sys.S =  sym('S', [sys.X_DIMS, sys.X_DIMS]);
sys.S = tril(sys.S,0) + tril(sys.S,-1).';
sys.a = sym('a', [sys.X_DIMS, 1]);
sys.b = sym('b', [sys.X_DIMS, 1]);
sys.err_lqr = sys.x.'*sys.S*sys.x - sys.x.'*S_joint*sys.x;
for ii=1:1:sys.X_DIMS
    sys.err_lqr = int(sys.err_lqr, sys.x(ii), [sys.a(ii), sys.b(ii)]);
end
sys.err_lqr_func = matlabFunction(simplify(sys.err_lqr), 'Vars', {sys.S, sys.a, sys.b});

sys.state_bounds = [-0.5, 0.5; -1, 1; 2*pi/3, 4*pi/3; -1, 1];
sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1));

err_lqr = zeros(1, size(u_x,1));
for dd=1:1:size(u_x,1)
    p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), [sys.U_DIMS, 2]);
    s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), [sys.U_DIMS, sys.X_DIMS]);
    err_lqr(dd) = computeLQRMeasure(sys, p, s);
end