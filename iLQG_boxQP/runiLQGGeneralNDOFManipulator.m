clear;
close all;
clc;

%%

restoredefaultpath;
system_name = 'manipulator2dof';
addpath('general');
addpath(strcat('new_systems/', system_name));
load(strcat('new_systems/', system_name, '/sys.mat'), 'sys');

sys.m = [5; 1]; % kg
sys.l = [0.24; 0.9]; % m
sys.g = 9.81; % m/s^2
sys.T = 4;
sys.dt = 0.001;
sys.X_DIMS = 2*sys.n; % [thi, ... dthi, ...]
sys.U_DIMS = sys.n;   % [taui] 

sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(1) = pi;
sys.goal = sys.l_point;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.fxu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
% sys.lims = 15*[-ones(sys.n,1), ones(sys.n,1)]; % action limits
sys.lims = 15*[-1, 1;
               -1, 1]; % action limits

sys.gamma_ = 0.997;
sys.lambda_ = (1 - sys.gamma_)/sys.dt;
sys.Q = diag([25*ones(1, sys.n), 0.02*ones(1, sys.n)]);
sys.R = diag(0.001*ones(1, sys.n));
sys.full_DDP = false;

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
% Op.Alpha = [1];

theta1_starts = [2*pi/3, 2*pi/3, 4*pi/3, 4*pi/3;
                 -1,     1,      -1,     1];
theta2_starts = [-pi/3, -pi/3, pi/3, pi/3;
                 -1,    1,     -1,   1];

starts = zeros(sys.X_DIMS, size(theta1_starts, 2) ...
                           *size(theta2_starts, 2));
count = 0;
for th1=1:1:size(theta1_starts,2)
    for th2=1:1:size(theta2_starts, 2)
        count = count + 1;
        starts([1;3], count) = theta1_starts(:, th1);
        starts([2;4], count) = theta2_starts(:, th2);
    end
end
starts = starts(:, [1, end]);

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
sys.name = 'manipulator2dofwu0';
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
sys.name = 'manipulator2dof';
p_joint = [0, 1;
           0, 1];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjoint, Ujoint, cjoint] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);

sys.u0init = true;
sys.name = 'manipulator2dofwu0';
p_joint = [0, 1;
           0, 1];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjointwu0, Ujointwu0, cjointwu0] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);
