clear;
close all;
clc;
%%

restoredefaultpath;
system_name = 'cartpole';
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/',system_name,'System.mat'));

assert(isfield(sys, 'name') && strcmp(sys.name, system_name), 'Check loaded system!');

sys.X_DIMS = 4;
sys.U_DIMS = 2;
sys.mc = 5;
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81; % m/s^2
sys.Q = diag([25, 0.02, 25, 0.02]);
sys.R = diag([0.001, 0.001]);
sys.gamma_ = 0.997;
sys.dt = 0.001;
sys.limits = [-1.5, 1.5;
              -3, 3;
              0, 2*pi;
              -3, 3];
sys.limits_t = [0, 4];
sys.lims = [-9, 9;
            -9, 9];

sys.l_point = [0; 0; pi; 0];
% Trajectory to track
traj_points = 5;
traj_start = [-1; 0; 0; 0];
traj_end = sys.l_point;
sys.goal = zeros(sys.X_DIMS, traj_points);
for xx=1:1:sys.X_DIMS
    sys.goal(xx, :) = linspace(traj_start(xx), traj_end(xx), traj_points);
end
sys.u0 = [0; 0];

Op.num_points = [31, 31, 31, 31];
Op.num_action_samples = [12, 12];
Op.max_iter = 1000;
Op.max_policy_iter = 100;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1))*2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1))/12;
Op.gtol = 0.000002;
Op.save_dir = 'data';
Op.reuse_policy = false;

% u_x = [];
% p = [0, 1;1, 1];
% s = [1, 0, 0, 0;
%      0, 1, 1, 1];
% 
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];
% 
% policies = cell(size(u_x,1), 1);
% value = cell(size(u_x,1), 1);
% info = cell(size(u_x,1), 1);
% 
% for dd=1:1:1
%     
%     p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
%     s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
%     
%     disp(sprintf('Decomposition %d/%d', dd, size(u_x,1)));
%     sys.decomposition_id = dd;
%     [policies{dd,1}, value{dd,1}, info{dd,1}] = trajdp_decomposition(sys, Op, p, s);
%     
% end

p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;

disp('Joint');
[policies_joint, value_joint, info_joint] = trajdp_decomposition(sys, Op, p_joint, s_joint);
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
