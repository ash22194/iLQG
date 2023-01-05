clear;
close all;
clc;
num_gpus = gpuDeviceCount();
gpu_id = 0;
max_avail_memory = 0;
for gg=1:1:num_gpus
    g = gpuDevice(gg);
    if (g.AvailableMemory > max_avail_memory)
        gpu_id = gg;
        max_avail_memory = g.AvailableMemory;
    end
end
g = gpuDevice(gpu_id);
reset(g);
fprintf('Using GPU : %d\n', gpu_id);

%%

restoredefaultpath;
system_name = 'cartpole';
addpath('../iLQG_boxQP');
addpath('../iLQG_boxQP/general');
addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
load(strcat('data/',system_name,'System.mat'));
mexcuda(strcat('systems/', system_name, '/dyn_mex_finite.cu'), '-R2018a', '-output', strcat('systems/', system_name, '/dyn_mex_finite'));

assert(isfield(sys, 'name') && strcmp(sys.name, system_name), 'Check loaded system!');

sys.X_DIMS = 4;
sys.U_DIMS = 2;
sys.mc = 5;
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81; % m/s^2
sys.Q = diag([25, 0.02, 25, 0.02]);
sys.Qf = 100 * sys.Q;
sys.R = diag([0.001, 0.001]);
sys.gamma_ = 0.997;
sys.dt = 0.001;

sys.limits_t = [0, 4];
% sys.limits_t = [0, 4];
sys.lims = [-9, 9;
            -9, 9];

sys.l_point = [0; 0; pi; 0];

sys.goal = sys.l_point;
% Trajectory to track
traj_points = 60;
traj_start = [-0.5; 0; pi/2; 0];
% traj_points = 80;
% traj_start = [-1; 0; 0; 0];
trajectory = cell(4,1);

Op_traj.lims  = sys.lims;
Op_traj.maxIter = 500;
Op_traj.gamma_ = sys.gamma_;
Op_traj.print = 0;
Op_traj.Alpha = [0.1];

sys.T = sys.limits_t(2);
sys.dt = sys.T / traj_points;
sys.X_DIMS_FREE = linspace(1, sys.X_DIMS, sys.X_DIMS)';
sys.X_DIMS_FIXED = [];
sys.U_DIMS_FREE = linspace(1, sys.U_DIMS, sys.U_DIMS)';
sys.U_DIMS_FIXED = [];
sys.U_DIMS_CONTROLLED = [];
sys.u0 = zeros(sys.U_DIMS, 1);
sys.u0init = true;
sys.full_DDP = false;
sys.cxmod = zeros(sys.X_DIMS, 1);
[trajectory{4}, ~, trajectory{1}, trajectory{2}, trajectory{3}, ...
 ~, ~, ~, ~, ~, ~, ~] = get_ilqg_trajectory_multitraj_parallel(sys, Op_traj, traj_start, [], ...
                                                               cell(0,6), cell(0,6));

restoredefaultpath;
addpath(strcat('systems/', system_name));
addpath('systems');

% 'Resample' as per the trajectory length
sys.goal = zeros(sys.X_DIMS, traj_points);
sys.u0 = zeros(sys.U_DIMS, traj_points);
sys.goal(:,1) = trajectory{1}(:,1);
dt_traj = sys.dt;
sys.limits = inf * [-ones(sys.X_DIMS,1), ones(sys.X_DIMS,1)];
for tt=1:1:(traj_points-1)
    [~, closest_x] = min(vecnorm(sys.goal(:,tt) - trajectory{4}, 2, 1));
    sys.u0(:, tt) = trajectory{2}(:,closest_x) ...
                    + trajectory{3}(:,:,closest_x) * (sys.goal(:,tt) - trajectory{1}(:,closest_x));
    sys.goal(:, tt+1) = cell2mat(dyn_finite_rk4(sys, num2cell(sys.goal(:, tt)), num2cell(sys.u0(:, tt)), dt_traj));
end
sys.goal(:,end)
sys.limits = [-2, 1;
              -3, 3;
              0, 2*pi;
              -3, 5];
% sys.limits = [-1.5, 1.5;
%               -3, 3;
%               0, 2*pi;
%               -3, 3];
sys.dt = 1e-3;

% traj_start = sys.l_point;
% traj_end = sys.l_point;
% sys.goal = zeros(sys.X_DIMS, traj_points);
% for xx=1:1:sys.X_DIMS
%     sys.goal(xx, :) = linspace(traj_start(xx), traj_end(xx), traj_points);
% end
% sys.u0 = zeros(sys.U_DIMS, traj_points);

Op.num_points = [25, 25, 25, 25];
% Op.num_points = [21, 21, 21, 21];

Op.num_action_samples = [12, 12];
Op.max_iter = 1000;
Op.max_policy_iter = 100;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1))*2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1))/12;
Op.gtol = 0.000002;
Op.save_dir = 'data/trajectory_tracking';
Op.reuse_policy = true;

X_DIMENSIONS = linspace(1,4,4);
u_x = [];
Cu = cell(0, 2);

% F first
p = [2, 1;0, 1];
p = reshape(p, 1, 4);
for jj=1:1:4
    s1 = nchoosek(X_DIMENSIONS, jj);
    for ss=1:1:size(s1, 1)
        s = [zeros(1, 4); ones(1, 4)];
        s(1,s1(ss,:)) = 1;
        s(2,s1(ss,:)) = 0;
        u_x = [u_x; [p, reshape(s, 1, 8)]];
        Cu = cat(1, Cu, {s1(ss,:), 1});
    end
end

% T first
p = [0, 1;1, 1];
p = reshape(p, 1, 4);
for jj=1:1:4
    s1 = nchoosek(X_DIMENSIONS, jj);
    for ss=1:1:size(s1, 1)
        s = [ones(1, 4); zeros(1, 4)];
        s(1,s1(ss,:)) = 0;
        s(2,s1(ss,:)) = 1;
        u_x = [u_x; [p, reshape(s, 1, 8)]];
        Cu = cat(1, Cu, {s1(ss,:), 2});
    end
end

% Decoupled
p = [0, 1;0, 2];
p = reshape(p, 1, 4);
for jj=1:1:3
    s1 = nchoosek(X_DIMENSIONS, jj);
    for ss=1:1:size(s1, 1)
        s = [zeros(1, 4); ones(1, 4)];
        s(1,s1(ss,:)) = 1;
        s(2,s1(ss,:)) = 0;
        u_x = [u_x; [p, reshape(s, 1, 8)]];
        Cu = cat(1, Cu, {s1(ss,:), 3});
    end
end

% p = [0, 1;
%      1, 1];
% s = [1,1,0,0;
%      0,0,1,1];
% u_x = [reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.X_DIMS*sys.U_DIMS)];

% s = [0,1,0,0;
%      1,0,1,1];
% u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.X_DIMS*sys.U_DIMS)];

% s = [1,0,0,0;
%      0,1,1,1];
% u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.X_DIMS*sys.U_DIMS)];

% s = [0,0,0,0;
%      1,1,1,1];
% u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.X_DIMS*sys.U_DIMS)];

% p = [2, 1;
%      0, 1];
% s = [1,1,1,0;
%      0,0,0,1];
% u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.X_DIMS*sys.U_DIMS)];

% s = [1,1,0,0;
%      0,0,1,1];
% u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.X_DIMS*sys.U_DIMS)];

sys.name = 'cartpole_finite_horizon';
value = cell(size(u_x,1), 1);

for dd=1:1:size(u_x,1)
    p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    disp(sprintf('Decomposition %d/%d', dd, size(u_x,1)));
    sys.decomposition_id = dd;
    [~, value{dd,1}, ~] = trajdp_decomposition_gpu(sys, Op, p, s);
    value{dd,1} = gather(value{dd,1});
end

p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;

disp('Joint');
[~, value_joint, info_joint] = trajdp_decomposition_gpu(sys, Op, p_joint, s_joint);
value_joint = gather(value_joint);
info_joint.state_grid = cellfun(@(x) gather(x), info_joint.state_grid, 'UniformOutput', false);

grid_indices = cell(sys.X_DIMS, 1);
for xxi=1:1:(sys.X_DIMS)
    grid_indices{xxi} = linspace(1, Op.num_points(xxi), Op.num_points(xxi));
end
sys.state_bounds = [-0.5, 0.5; -1, 1; -pi/3, pi/3; -1, 1];
valid_range_final = true(size(info_joint.state_grid{1}));
err_dp = zeros(size(u_x,1), traj_points);
for tt=1:1:traj_points
    state_bounds = mat2cell(sys.state_bounds + sys.goal(:,tt), ones(sys.X_DIMS,1), 2);
    valid_range = cellfun(@(x,y) (x(grid_indices{:}, tt) >= y(1)) & (x(grid_indices{:}, tt) <= y(2)), ...
                          info_joint.state_grid(1:sys.X_DIMS), state_bounds, 'UniformOutput', false);
    for dim=1:1:(sys.X_DIMS)
        valid_range_final(grid_indices{:}, tt) = and(valid_range_final(grid_indices{:}, tt), ...
                                                     valid_range{dim, 1});
    end
    
    for dd=1:1:size(u_x,1)
        value_d = value{dd,1}(grid_indices{:}, tt);
        value_j = value_joint(grid_indices{:}, tt);
        err_dp(dd, tt) = mean(abs(value_d(valid_range_final(grid_indices{:}, tt)) ...
                                  - value_j(valid_range_final(grid_indices{:}, tt))), 'all');
    end
end

save(strcat(Op.save_dir, '/', system_name, '/summary.mat'), 'u_x', 'value', 'value_joint', 'info_joint', 'err_dp', 'sys', 'Op', '-v7.3', '-nocompression');
