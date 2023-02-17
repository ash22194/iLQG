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
Op.num_action_samples = [10,4,4,10];
Op.max_iter = 5000;
Op.max_policy_iter = 500;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1))*2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1))/12;
Op.gtol = 0.00000002*0;
Op.save_dir = 'data';
Op.save_every = 500;
Op.reuse_policy = true;

u_x = [];

% GA
% p = [3, 1;3, 1;0, 1;0, 2];
% s = [1, 1, 0, 0, 0, 1, 1, 1, 1, 0;
%      1, 1, 0, 0, 0, 1, 1, 1, 1, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 0, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% GA
% p = [2, 1;0, 1;0, 1;0, 2];
% s = [1, 1, 1, 0, 0, 0, 1, 0, 1, 0;
%      0, 0, 0, 0, 1, 1, 0, 1, 0, 0;
%      0, 0, 0, 0, 1, 1, 0, 1, 0, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

p = [2, 1;0, 1;0, 1;0, 2];
s = [1, 0, 0, 0, 0, 0, 1, 1, 0, 0;
     0, 1, 1, 0, 1, 1, 0, 0, 1, 0;
     0, 1, 1, 0, 1, 1, 0, 0, 1, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% MCTS
p = [0, 2;0, 2;0, 2;0, 1];
s = [1, 1, 1, 0, 1, 1, 1, 1, 1, 0;
     1, 1, 1, 0, 1, 1, 1, 1, 1, 0;
     1, 1, 1, 0, 1, 1, 1, 1, 1, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Random % Too memory intensive to evaluate
% p = [2, 1;0, 1;0, 1;2, 2];
% s = [1, 0, 0, 0, 0, 1, 1, 1, 0, 0;
%      0, 1, 1, 0, 0, 0, 0, 0, 1, 0;
%      0, 1, 1, 0, 0, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 1, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Pareto : F1 - (z, vz), F2 - (r, vr, vy, F1), F3 - (p, vx, vp, F1, F2), F4 - (ya, vya)
% p = [2, 2; 3, 1; 0, 3; 0, 4];
% s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
%      0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
%      0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
% u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Pareto
p = [3, 2; 0, 1; 2, 2; 0, 2];
s = [1, 0, 0, 0, 0, 0, 1, 0, 1, 0;
     0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
     0, 0, 1, 0, 1, 0, 0, 0, 0, 0;
     0, 0, 0, 1, 0, 0, 0, 0, 0, 1];
u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

% Self created
% p = [0, 1; 0, 2; 0, 3; 0, 4];
% s = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0;
%      0, 1, 0, 0, 0, 1, 0, 1, 0, 0;
%      0, 0, 1, 0, 1, 0, 0, 0, 1, 0;
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
    disp(strcat('Computation time: ', num2str(info{dd,1}.time_total)));
end
% return;
decomposition = load(strcat(Op.save_dir, '/', sys.name, '/decomp4/final.mat'));
policies = cat(1, policies, {decomposition.policies});
info = cat(1, info, {decomposition.info});
value = cat(1, value, {0});

% p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
% s_joint = ones(sys.U_DIMS, sys.X_DIMS);
% sys.decomposition_id = 0;

% disp('Joint');
% [policies_joint, value_joint, info_joint] = dp_decomposition_gpu(sys, Op, p_joint, s_joint);
% policies_joint = cellfun(@(x) gather(x), policies_joint, 'UniformOutput', false);
% value_joint = gather(value_joint);

state_bounds = [0.7000, 1.3000;
               -0.7854, 0.7854;
               -0.7854, 0.7854;
               -0.7854, 0.7854;
               -1.0000, 1.0000;
               -1.0000, 1.0000;
               -0.7500, 0.7500;
               -0.5000, 0.5000;
               -0.5000, 0.5000;
               -0.2500, 0.2500];
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

% Closed-loop trajectories
policy_interpolant = cell(size(policies, 1), 1);
policy_function = cell(size(policies, 1), 1);
for dd = 1:1:size(policies,1)
     policy_interpolant{dd, 1} = cell(sys.U_DIMS, 1);
     policy_function{dd, 1} = cell(sys.U_DIMS, 1);
     for uu=1:1:sys.U_DIMS
          uui = find(cellfun(@(x) any(x==uu), policies{dd, 1}(:,1)));
          uuii = find(policies{dd, 1}{uui, 1}==uu);
          policy_interpolant{dd, 1}{uu} = griddedInterpolant(policies{dd, 1}{uui, 4}.state_grid{:}, policies{dd, 1}{uui, 3}{uuii, 1});
          policy_function{dd, 1}{uu} = @(x) policy_interpolant{dd, 1}{uu}(x(policies{dd, 1}{uui, 2}));
     end
end
% policy_interpolant{size(u_x,1)+1, 1} = cellfun(@(x) griddedInterpolant(info_joint.state_grid{:}, x), ...
%                                                policies_joint{1,3}, 'UniformOutput', false);
% policy_function{size(u_x,1)+1, 1} = cellfun(@(x) @(y) x(y(policies_joint{1,2})), policy_interpolant{size(u_x,1)+1, 1}, 'UniformOutput', false);

% state_bounds = cell2mat(state_bounds);
time_span = [0, 12];
NUM_CTRL = round(time_span(2) / sys.dt);
num_starts = 100;
starts = sys.l_point + 0.9*repmat(state_bounds(:,2) - state_bounds(:,1), [1, num_starts]).*(rand(sys.X_DIMS, num_starts) - 0.5);
trajectories = cell(size(policy_function, 1), num_starts, 2);
costs = zeros(num_starts, size(policy_function, 1));
% load(strcat('data/', system_name, '/trajectories_Pareto_GA_MCTS_Heuristic.mat'), 'starts');
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

save(strcat('data/', system_name, '/trajectories_GA_MCTS_Pareto_Heuristic_100.mat'), 'starts', 'trajectories', 'time_span', 'sys', 'Op', 'costs', 'u_x', '-v7.3');

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
