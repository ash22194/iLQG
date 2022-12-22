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
system_name = 'biped2d';
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/',system_name,'System.mat'));
mexcuda(strcat('systems/', system_name, '/dyn_mex_finite.cu'), '-R2018a', '-output', strcat('systems/', system_name, '/dyn_mex_finite'));

assert(isfield(sys, 'name') && strcmp(sys.name, system_name), 'Check loaded system!');

sys.X_DIMS = 6;
sys.U_DIMS = 4;
sys.m = 72;
sys.I = 3;
sys.d = 0.2;
sys.df = 0.5;
sys.l0 = 1.15;
sys.g = 9.81; % m/s^2
sys.Q = diag([350, 700, 1.5, 1.5, 500, 5]);
sys.R = diag([0.000002*ones(1, 2),  0.00002*ones(1, 2)])/2;
sys.gamma_ = 0.999;
sys.dt = 0.001;
sys.limits = [sys.l0 - 0.3, sys.l0 + 0.1;
              pi/2, pi/2 + 0.6;
              -0.3, 0.5;
              -0.5, 1;
              -pi/8, pi/8;
              -2, 2];
sys.lims = 2*[0, 1.5*sys.m*sys.g;
              0, 1.5*sys.m*sys.g;
              -0.125*sys.m*sys.g, 0.125*sys.m*sys.g;
              -0.125*sys.m*sys.g, 0.125*sys.m*sys.g];
lg = 0.96;
alpha1g = pi/2 + asin(sys.df/2/lg);
l2g = sqrt((sys.df + lg*cos(alpha1g))^2 + (lg*sin(alpha1g))^2);
alpha2g = acos((sys.df + lg*cos(alpha1g))/l2g);
sys.l_point = [lg; alpha1g; 0; 0; 0; 0];
sys.goal = sys.l_point;
sys.u0 = [sys.m*sys.g*cos(alpha2g)/sin(alpha1g - alpha2g);
          -sys.m*sys.g*cos(alpha1g)/sin(alpha1g - alpha2g);
          0;
          0];

Op.num_points = [13,13,14,19,14,21];
Op.num_action_samples = [5,5,4,4];
Op.max_iter = 2000;
Op.max_policy_iter = 100;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1))*2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1))/12;
Op.gtol = 0.000002;
Op.save_dir = 'data';
Op.reuse_policy = true;

% u_x = [];
% % F - COM, T - All
% p = [3,1; 3,1; 0,1; 0,1];
% s = [ones(2,4), zeros(2,2);
%      zeros(2,4), ones(2,2)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % T - Torso, F - All
% p = [0,1; 0,1; 1,1; 1,1];
% s = [ones(2,4), zeros(2,2);
%      zeros(2,4), ones(2,2)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % T - All, F - All
% p = [0,1; 0,1; 1,1; 1,1];
% s = [zeros(2,6);
%      ones(2,6)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % F - All, T - All
% p = [3,1; 3,1; 0,1; 0,1];
% s = [ones(2,6);
%      zeros(2,6)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % T - Torso, F - COM
% p = [0,1; 0,1; 0,2; 0,2];
% s = [ones(2,4), zeros(2,2);
%      zeros(2,4), ones(2,2)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % F - Torso, T - All
% p = [3,1; 3,1; 0,1; 0,1];
% s = [zeros(2,4), ones(2,2);
%      ones(2,4), zeros(2,2)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % T - COM, F - All
% p = [0,1; 0,1; 1,1; 1,1];
% s = [zeros(2,4), ones(2,2);
%      ones(2,4), zeros(2,2)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % T - COM, F - Torso
% p = [0,1; 0,1; 0,2; 0,2];
% s = [zeros(2,4), ones(2,2);
%      ones(2,4), zeros(2,2)];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

u_x = [];
% GA
% 1
p = [0, 1;
     0, 2;
     0, 3;
     1, 1];
s = [0, 0, 0, 1, 0, 0;
     0, 1, 1, 0, 0, 0;
     0, 0, 0, 0, 1, 1;
     1, 0, 0, 0, 0, 0];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% MCTS
p = [0, 1;
     0, 2;
     1, 1;
     0, 3];
s = [0, 0, 0, 1, 0, 0;
     0, 1, 1, 0, 0, 0;
     1, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% p = [0, 1;
%      0, 2;
%      0, 3;
%      0, 3];
% s = [1, 0, 0, 1, 0, 0;
%      0, 1, 1, 0, 0, 0;
%      0, 0, 0, 0, 1, 1;
%      0, 0, 0, 0, 1, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% p = [0, 1;
%      0, 2;
%      1, 1;
%      0, 3];
% s = [1, 0, 0, 0, 0, 0;
%      0, 1, 1, 0, 0, 0;
%      0, 0, 0, 1, 0, 0;
%      0, 0, 0, 0, 1, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Random
% p = [0, 1;
%      0, 2;
%      0, 3;
%      2, 1];
% s = [0, 1, 1, 0, 0, 0;
%      1, 0, 0, 0, 0, 0;
%      0, 0, 0, 0, 1, 1;
%      0, 0, 0, 1, 0, 0];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

p = [0, 1;
     0, 2;
     0, 3;
     2, 1];
s = [1, 0, 0, 1, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 1, 1;
     0, 1, 0, 0, 0, 0];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% p = [0, 1;
%      3, 1;
%      0, 2;
%      0, 3];
% s = [1, 0, 0, 1, 0, 0;
%      0, 1, 1, 0, 0, 0;
%      0, 0, 0, 0, 0, 0;
%      0, 0, 0, 0, 1, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% p = [0, 1;
%      0, 2;
%      0, 3;
%      0, 3];
% s = [1, 0, 0, 1, 0, 0;
%      0, 1, 1, 0, 0, 0;
%      0, 0, 0, 0, 1, 1;
%      0, 0, 0, 0, 1, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% % Pareto
% p = [0, 1;
%      0, 1;
%      0, 3;
%      0, 3];
% s = [1, 1, 1, 1, 0, 0;
%      1, 1, 1, 1, 0, 0;
%      0, 0, 0, 0, 1, 1;
%      0, 0, 0, 0, 1, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% F - COM, T - All
p = [3,1; 3,1; 0,1; 0,1];
s = [ones(2,4), zeros(2,2);
     zeros(2,4), ones(2,2)];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

policies = cell(size(u_x,1), 1);
value = cell(size(u_x,1), 1);
info = cell(size(u_x,1), 1);

for dd=1:1:size(u_x, 1)
    p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    sys.name = 'biped2d';
    disp(sprintf('Decomposition %d/%d', dd, size(u_x,1)));
    sys.decomposition_id = dd;
    [policies{dd,1}, value{dd,1}, info{dd,1}] = dp_decomposition_gpu(sys, Op, p, s);
    policies{dd,1} = cellfun(@(x) gather(x), policies{dd,1}, 'UniformOutput', false);
    value{dd,1} = gather(value{dd,1});
    info{dd,1}.state_grid = cellfun(@(x) gather(x), info{dd,1}.state_grid, 'UniformOutput', false);
end

p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;

disp('Joint');
sys.name = 'biped2d/decomp_joint';
[policies_joint, value_joint, info_joint] = dp_decomposition_gpu(sys, Op, p_joint, s_joint);
policies_joint = cellfun(@(x) gather(x), policies_joint, 'UniformOutput', false);
value_joint = gather(value_joint);
info_joint.state_grid = cellfun(@(x) gather(x), info_joint.state_grid, 'UniformOutput', false);

state_bounds = [0.95, 1;
                0.3, 0.4;
                -0.1, 0.1;
                -0.3, 0.3;
                -0.2, 0.2;
                -0.2, 0.2];
state_bounds(2,:) = state_bounds(2,:) + pi/2;
state_bounds = mat2cell(state_bounds, ones(sys.X_DIMS,1), 2);

valid_range = cellfun(@(x,y) (x>y(1)) & (x<y(2)), info_joint.state_grid, state_bounds, 'UniformOutput', false);
valid_range_final = true(size(info_joint.state_grid{1}));
for dim=1:1:(sys.X_DIMS)
    valid_range_final = and(valid_range_final, valid_range{dim,1});
end
valid_range = valid_range_final;

err_dp = cellfun(@(x) mean(x(valid_range) - value_joint(valid_range), 'all'), value);

% Trajectory based estimate
policy_interpolant = cell(size(u_x, 1)+1, 1);
policy_function = cell(size(u_x, 1)+1, 1);
for dd = 1:1:size(u_x,1)
     policy_interpolant{dd, 1} = cell(sys.U_DIMS, 1);
     policy_function{dd, 1} = cell(sys.U_DIMS, 1);
     for uu=1:1:sys.U_DIMS
          uui = find(cellfun(@(x) any(x==uu), policies{dd, 1}(:,1)));
          uuii = find(policies{dd, 1}{uui, 1}==uu);
          policy_interpolant{dd, 1}{uu} = griddedInterpolant(policies{dd, 1}{uui, 4}.state_grid{:}, policies{dd, 1}{uui, 3}{uuii, 1});
          policy_function{dd, 1}{uu} = @(x) policy_interpolant{dd, 1}{uu}(x(policies{dd, 1}{uui, 2}));
     end
end
policy_interpolant{size(u_x,1)+1, 1} = cellfun(@(x) griddedInterpolant(info_joint.state_grid{:}, x), ...
                                               policies_joint{1,3}, 'UniformOutput', false);
policy_function{size(u_x,1)+1, 1} = cellfun(@(x) @(y) x(y(policies_joint{1,2})), policy_interpolant{size(u_x,1)+1, 1}, 'UniformOutput', false);

state_bounds = cell2mat(state_bounds);
time_span = [0, 12];
NUM_CTRL = round(time_span(2) / sys.dt);
num_starts = 50;
starts = state_bounds(:,1) + repmat(state_bounds(:,2) - state_bounds(:,1), [1, num_starts]).*rand(sys.X_DIMS, num_starts);
trajectories = cell(size(policy_function, 1), num_starts, 2);
costs = zeros(num_starts, size(policy_function, 1));
sys.X_DIMS_FREE = 1:1:sys.X_DIMS;
for dd=1:1:size(policy_function, 1)
     disp(strcat('Decomposition :', num2str(dd)));
     for nn = 1:1:num_starts
          start = starts(:, nn);
          trajectories{dd, nn, 1} = zeros(sys.X_DIMS, NUM_CTRL+1);
          trajectories{dd, nn, 1}(:,1) = start;
          trajectories{dd, nn, 2} = sys.u0 * ones(1, NUM_CTRL+1);
          discount = 1;
          for tt = 1:1:NUM_CTRL
               x = trajectories{dd, nn, 1}(:, tt);
               [xnext, u] = dyn_closed_loop_rk4(sys, x, policy_function{dd, 1}, sys.dt);
               costs(nn, dd) = costs(nn, dd) + sys.dt*discount*((x - sys.goal)'*sys.Q*(x - sys.goal) ...
                                                                + (u - sys.u0)'*sys.R*(u - sys.u0));
               discount = discount * sys.gamma_;
               trajectories{dd, nn, 1}(:, tt+1) = xnext;
               trajectories{dd, nn, 2}(:, tt) = u;
          end
          disp(strcat('Trajectory ', num2str(nn), ', End :', sprintf('%f ', trajectories{dd, nn, 1}(:,end))));
     end
end
err_dp_traj = mean(costs(:,1:(size(costs,2) - 1)) - costs(:, end), 1);

% save(strcat(Op.save_dir, '/', system_name, '/summary.mat'), 'u_x', 'policies', 'value', 'info', 'policies_joint', 'value_joint', 'info_joint', 'err_dp', 'sys', 'Op', '-v7.3', '-nocompression');
% save(strcat(Op.save_dir, '/', system_name, '/summary_GA_MCTS_Random.mat'), 'u_x', 'policies', 'value', 'info', 'sys', 'Op', '-v7.3', '-nocompression');

function [xnext, u] = dyn_closed_loop_rk4(sys, x, policies, dt)

     x_ = max(sys.limits(:,1), min(sys.limits(:,2), x(:)));
     x_ = num2cell(x_);
     x = num2cell(x);

     u = num2cell(cellfun(@(y) y(x_), policies));
     k1 = dyn(sys, x, u);
     q = num2cell(cellfun(@(y,z) y + z*dt*0.5, x, k1));

     k2 = dyn(sys, q, u);
     q = num2cell(cellfun(@(y,z) y + z*dt*0.5, x, k2));

     k3 = dyn(sys, q, u);
     q = num2cell(cellfun(@(y,z) y + z*dt, x, k3));

     k4 = dyn(sys, q, u);

     xnext = cell2mat(x) + dt/6*(cell2mat(k1) + 2*cell2mat(k2) + 2*cell2mat(k3) + cell2mat(k4));
     u = cell2mat(u);
end