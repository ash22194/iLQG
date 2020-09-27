clear;
close all;
clc;

%%

restoredefaultpath;
system_name = 'cartpole';
addpath('general');
addpath(strcat('new_systems/', system_name));
load(strcat('new_systems/', system_name, '/sys.mat'), 'sys');

sys.mc = 5;
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81;
sys.T = 4;
sys.dt = 0.001;
sys.X_DIMS = 4; % x, x_dot, th, th_dot
sys.U_DIMS = 2; % F, tau 

sys.l_point = [0;0;pi;0];
sys.goal = sys.l_point;
sys.u0 = [0;0];
sys.fxu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
sys.lims = 9*[-ones(2,1), ones(2,1)];

sys.gamma_ = 0.997;
sys.lambda_ = (1 - sys.gamma_)/sys.dt;
sys.Q = diag([25, 0.02, 25, 0.02]);
sys.R = 0.001*eye(2);
sys.full_DDP = false;

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
% Op.Alpha = [1];

cart_starts = [-0.4, -0.2, 0.2, 0.4;
                0,     0,  0,   0];
pole_starts = [2*pi/3, 3*pi/4, 5*pi/4, 4*pi/3;
               0,      0,      0,      0];
starts = zeros(sys.X_DIMS, size(cart_starts, 2) ...
                           *size(pole_starts, 2));
count = 0;
for carts=1:1:size(cart_starts, 2)
    for poles=1:1:size(pole_starts, 2)
        count = count + 1;
        starts(1:2, count) = cart_starts(:, carts);
        starts(3:4, count) = pole_starts(:, poles);
    end
end
starts = starts(:,[1, end]);

%% Test Decompositions

decompositions_file = strcat('../GA/data/', system_name, 'ParetoFront.mat');
load(decompositions_file, 'u_x', 'u_err_lqr');
[~, u_min_id] = min(u_err_lqr(:,1));
u_x = u_x(u_min_id, :);

sys.u0init = false;
Xd = zeros(sys.X_DIMS, round(sys.T / sys.dt)+1, size(starts, 2), size(u_x, 1));
Ud = zeros(sys.U_DIMS, round(sys.T / sys.dt), size(starts, 2), size(u_x, 1));
cd = zeros(1, size(starts, 2), size(u_x, 1));

for d=1:1:size(u_x, 1)
    
    fprintf('Decomposition : %d / %d\n', d, size(u_x, 1));
    sys.decomposition_id = d;
    p = reshape(u_x(d, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(d, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    [Xd(:,:,:,d), Ud(:,:,:,d), cd(:,:,d)] = ilqg_decomposition(sys, Op, p, s, starts);
    
end

sys.u0init = true;
sys.name = 'cartpolewu0';
Xdwu0 = zeros(sys.X_DIMS, round(sys.T / sys.dt)+1, size(starts, 2), size(u_x, 1));
Udwu0 = zeros(sys.U_DIMS, round(sys.T / sys.dt), size(starts, 2), size(u_x, 1));
cdwu0 = zeros(1, size(starts, 2), size(u_x, 1));

for d=1:1:size(u_x, 1)
    
    fprintf('Decomposition : %d / %d\n', d, size(u_x, 1));
    sys.decomposition_id = d;
    p = reshape(u_x(d, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(d, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    [Xdwu0(:,:,:,d), Udwu0(:,:,:,d), cdwu0(:,:,d)] = ilqg_decomposition(sys, Op, p, s, starts);
    
end

%% Test Joint Optimization

sys.u0init = false;
sys.name = 'cartpole';
p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjoint, Ujoint, cjoint] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);

sys.u0init = true;
sys.name = 'cartpolewu0';
p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjointwu0, Ujointwu0, cjointwu0] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);


