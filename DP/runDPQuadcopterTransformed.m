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
system_name = 'quadcopter_transformed';
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/',system_name,'System.mat'));
mexcuda(strcat('systems/', system_name, '/dyn_mex_finite.cu'), '-R2018a', '-output', strcat('systems/', system_name, '/dyn_mex_finite')); 

% assert(isfield(sys, 'name') && strcmp(sys.name, system_name), 'Check loaded system!');

sys.X_DIMS = 10; % z, ro, pi, ya, vx, vy, vz, vr0, vpi, vya
sys.U_DIMS = 4;
sys.m = 0.5;
sys.I = diag([4.86*1e-3; 4.86*1e-3; 8.8*1e-3]);
sys.l = 0.225;
sys.g = 9.81;
sys.bk = 1.14*1e-7/(2.98*1e-6); % tau/f
sys.T = 6;
sys.dt = 0.00025;
% sys.Q = diag([10, 2, 2, 0.25, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01]);
% sys.R = 0.0025*eye(4);
% sys.Q = diag([1, 2, 2, 0.25, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01]); % Original
% sys.Q = diag([5, 2, 2, 1, 0.025, 0.025, 0.01, 0.01, 0.01, 0.01]); % Decent
% sys.Q = diag([5, 2.5, 2.5, 5, 0.05, 0.05, 0.01, 0.01, 0.01, 0.01]); % XY Drift
% sys.Q = diag([0.2, 5, 5, 0.2, 0.1, 0.1, 0.001, 0.1, 0.1, 0.001]);
% sys.Q = diag([5, 0.01, 0.01, 5, 0.1, 0.1, 0.05, 0.1, 0.1, 0.05]);
sys.Q = diag([5, 0.001, 0.001, 5, 0.5, 0.5, 0.05, 0.075, 0.075, 0.05]);
sys.R = diag([0.002, 0.01, 0.01, 0.004]);
sys.gamma_ = 0.99975;
% sys.limits = [0.4, 1.6;
%               -0.8, 0.8;
%               -0.8, 0.8;
%               -pi, pi;
%               -8, 8;
%               -8, 8;
%               -6, 6;
%               -11, 11;
%               -11, 11;
%               -6, 6];
% sys.limits = [0.4, 1.6;
%               -4*pi/9, 4*pi/9;
%               -4*pi/9, 4*pi/9;
%               -pi, pi;
%               -4, 4;
%               -4, 4;
%               -3, 3;
%               -6, 6;
%               -6, 6;
%               -4, 4];
sys.limits = [0.5, 1.5;
              -0.7, 0.7;
              -0.7, 0.7;
              -pi, pi;
              -2, 2;
              -2, 2;
              -1.5, 1.5;
              -6, 6;
              -6, 6;
              -2.5, 2.5];
sys.lims = [0, 2*sys.m*sys.g;
            -0.25*sys.m*sys.g, 0.25*sys.m*sys.g;
            -0.25*sys.m*sys.g, 0.25*sys.m*sys.g;
            -0.125*sys.m*sys.g, 0.125*sys.m*sys.g];
sys.l_point = [1; zeros(9,1)];
sys.goal = sys.l_point;
sys.u0 = [sys.m*sys.g; 0; 0; 0];

% Op.num_points = [37,37,37,75,81,81,81,111,111,81];
% Op.num_points = [17,17,17,35,25,25,25,51,51,35];
Op.num_points = [7,7,7,35,7,7,7,11,11,35];
% Op.num_points = 24*ones(1,sys.X_DIMS);
Op.num_action_samples = [10,10,10,10];
Op.max_iter = 5000;
Op.max_policy_iter = 800;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1))*2e-6/5;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1))/12/5;
Op.gtol = 0.00000002*0;
Op.save_dir = 'data';
Op.reuse_policy = false;

u_x = [];
% Pareto : F1 - (z, vz), F2 - (r, vr, vy, F1), F3 - (p, vx, vp, F1, F2), F4 - (ya, vya)
p = [2, 2; 3, 1; 0, 3; 0, 4];
s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
     0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
     0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Pareto : F1 - (z, vz), F2 - (r, vr, vy, F1), F3 - (p, vx, vp), F4 - (ya, vya)
% p = [2, 1; 0, 1; 0, 2; 0, 4];
% s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
%      0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Pareto : F1 - (z, vz), F2 - (r, vr, vy), F3 - (p, vx, vp), F4 - (ya, vya)
% p = [0, 1; 0, 2; 0, 3; 0, 4];
% s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
%      0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Self created
% p = [0, 1; 1, 1; 1, 2; 0, 2];
% s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
%      0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% GA
p = [3, 1;3, 1;0, 1;0, 2];
s = [1, 1, 0, 0, 0, 1, 1, 1, 1, 0;
     1, 1, 0, 0, 0, 1, 1, 1, 1, 0;
     0, 0, 1, 0, 1, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% GA
% p = [3, 1;3, 1;0, 1;0, 2];
% s = [1, 1, 0, 0, 0, 1, 1, 1, 0, 0;
%      1, 1, 0, 0, 0, 1, 1, 1, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% MCTS
p = [0, 2;0, 2;0, 2;0, 1];
s = [1, 1, 1, 0, 1, 1, 1, 1, 1, 0;
     1, 1, 1, 0, 1, 1, 1, 1, 1, 0;
     1, 1, 1, 0, 1, 1, 1, 1, 1, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % Self created
p = [0, 1; 0, 2; 0, 3; 0, 4];
s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
     0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
     0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % Self created
% p = [2, 1;0, 1;0, 1;0, 2];
% s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
%      0, 1, 1, 0, 1, 1, 0, 1, 1, 0;
%      0, 1, 1, 0, 1, 1, 0, 1, 1, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % Pareto
% p = [3, 2;0, 2;2, 1;0, 1];
% s = [1, 0, 1, 0, 1, 0, 1, 1, 0, 0;
%      0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
%      0, 1, 0, 0, 0, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

policies = cell(size(u_x,1), 1);
value = cell(size(u_x,1), 1);
info = cell(size(u_x,1), 1);
for dd=1:1:size(u_x,1)
    p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    sys.name = 'quadcopter_transformed';
    disp(sprintf('Decomposition %d/%d', dd, size(u_x,1)));
    sys.decomposition_id = dd;
    [policies{dd,1}, value{dd,1}, info{dd,1}] = dp_decomposition_gpu(sys, Op, p, s);
    policies{dd,1} = cellfun(@(x) gather(x), policies{dd,1}, 'UniformOutput', false);
    value{dd,1} = gather(value{dd,1});
    info{dd,1}.state_grid = cellfun(@(x) gather(x), info{dd,1}.state_grid, 'UniformOutput', false);
end

% p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
% s_joint = ones(sys.U_DIMS, sys.X_DIMS);
% sys.decomposition_id = 0;

% disp('Joint');
% [policies_joint, value_joint, info_joint] = dp_decomposition_gpu(sys, Op, p_joint, s_joint);
% policies_joint = cellfun(@(x) gather(x), policies_joint, 'UniformOutput', false);
% value_joint = gather(value_joint);
% info_joint.state_grid = cellfun(@(x) gather(x), info_joint.state_grid, 'UniformOutput', false);

% state_bounds = [0.95, 1;
%                 0.3, 0.4;
%                 -0.1, 0.1;
%                 -0.3, 0.3;
%                 -0.2, 0.2;
%                 -0.2, 0.2];
% state_bounds(2,:) = state_bounds(2,:) + pi/2;
% state_bounds = mat2cell(state_bounds, ones(sys.X_DIMS,1), 2);

% valid_range = cellfun(@(x,y) (x>y(1)) & (x<y(2)), info_joint.state_grid, state_bounds, 'UniformOutput', false);
% valid_range_final = true(size(info_joint.state_grid{1}));
% for dim=1:1:(sys.X_DIMS)
%     valid_range_final = and(valid_range_final, valid_range{dim,1});
% end
% valid_range = valid_range_final;

% err_dp = zeros(1, size(u_x,1));
% for dd=1:1:size(u_x,1)
%     err_dp(dd) = mean(abs(value{dd,1}(valid_range) - value_joint(valid_range)), 'all');
% end

% save(strcat(Op.save_dir, '/', system_name, '/summary.mat'), 'u_x', 'policies', 'value', 'info', 'policies_joint', 'value_joint', 'info_joint', 'err_dp', 'sys', 'Op', '-v7.3', '-nocompression');
% save(strcat(Op.save_dir, '/', system_name, '/summary_GA_MCTS_Random.mat'), 'u_x', 'policies', 'value', 'info', 'sys', 'Op', '-v7.3', '-nocompression');

