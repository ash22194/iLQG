clear;
close all;
% clc;
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
n = 4;
system_name = sprintf('manipulator%ddof', n);
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/',system_name,'System.mat'));
mexcuda(strcat('systems/', system_name, '/dyn_mex_finite.cu'), '-R2018a', '-output', strcat('systems/', system_name, '/dyn_mex_finite')); 

assert(isfield(sys, 'name') && strcmp(sys.name, system_name), 'Check loaded system!');

sys.X_DIMS = 2*sys.n; % [thi, ... dthi, ...]
sys.U_DIMS = sys.n;   % [taui]
if (n==2)
    sys.m = [2.5; 0.5]/2; % kg
    sys.l = [0.5; 0.25]/2; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8, 8, 0.6, 0.6])/5;
    sys.R = diag(0.003*(Izz(1)./Izz).^2);
    sys.limits = [0, 2*pi; repmat([-pi, pi], n-1, 1); repmat([-3, 3], n, 1)];
    sys.lims = 5*[-Izz/Izz(1), Izz/Izz(1)]; % action limits
    
    Op.num_points = 31 * ones(1, sys.X_DIMS);
    Op.num_action_samples = [15, 5];
    
    % Define decompositions to test
    u_x = [];
    
    % Cascaded
    p = [0, 1;1, 1];
    s = [1, 0, 1, 0;0, 1, 0, 1];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
    p = [0, 1;1, 1];
    s = [0, 1, 0, 1;1, 0, 1, 0];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
    p = [0, 1;1, 1];
    s = [0, 0, 0, 0;1, 1, 1, 1];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
    p = [2, 1;0, 1];
    s = [1, 0, 1, 0;0, 1, 0, 1];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
    p = [2, 1;0, 1];
    s = [0, 1, 0, 1;1, 0, 1, 0];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
    p = [2, 1;0, 1];
    s = [1, 1, 1, 1;0, 0, 0, 0];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
    % Decoupled
    p = [0, 1;0, 2];
    s = [1, 0, 1, 0;0, 1, 0, 1];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
    p = [0, 1;0, 2];
    s = [0, 1, 0, 1;1, 0, 1, 0];
    u_x = [u_x; reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
        
elseif (n==3)
    sys.m = [2.5; 0.5; 0.1] * 1.1; % kg
    sys.l = [0.5; 0.25; 0.125]; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8*ones(1,3), 0.6*ones(1,3)])/5;
    sys.R = diag(0.004*(Izz(1)./Izz));
    sys.limits = [0, 2*pi; repmat([-pi, pi], n-1, 1); repmat([-3, 3], n, 1)];
    sys.lims = [-16, 16; -7.5, 7.5; -1, 1]; % action limits
    
    Op.num_points = [17, 17, 17, 13, 13, 13];
    Op.num_action_samples = [8, 3, 2];
    
    % Define decompositions to test
%     load('data/manipulator3dof/manipulator3dof_paretofront.mat');
%     u_x = u_xp;
    p = [3, 1;
         3, 1;
         0, 1];
    s = [1, 1, 0, 1, 1, 0;
         1, 1, 0, 1, 1, 0;
         0, 0, 1, 0, 0, 1];
    u_x = [reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
elseif (n==4)
    sys.m = [5.4; 1.8; 0.6; 0.2]; % kg
    sys.l = [0.2; 0.5; 0.25; 0.125]; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8*ones(1,4), 0.2*ones(1,4)])/2;
    sys.R = diag([0.002; 0.004; 0.024; 0.1440]);
    sys.limits = [pi/2, 3*pi/2; repmat([-pi/2, pi/2], n-2, 1); [-3*pi/4, 3*pi/4]; repmat([-1.5, 1.5], n-1, 1); [-1.7, 1.7]];
    sys.lims = [-24, 24; -15, 15; -7.5, 7.5; -1, 1]; % action limits
    
    Op.num_points = [17,17,17,17,13,13,13,13];
    Op.num_action_samples = [6,4,3,3]*2;
    
    % Define decompositions to test
%     load('data/manipulator4dof/manipulator4dof_paretofront.mat');
%     u_x = u_xp;
    u_x = [];
    % GA
    p = [0, 1;
         0, 2;
         0, 3;
         0, 4];
    s = [eye(4), eye(4)];
    u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];
    
    % MCTS
%     p = [0, 1;
%          0, 2;
%          0, 3;
%          0, 3];
%     s = [1, 0, 0, 1, 1, 0, 0, 0;
%          0, 1, 0, 0, 0, 1, 0, 0;
%          0, 0, 1, 0, 0, 0, 1, 1;
%          0, 0, 1, 0, 0, 0, 1, 1];
%     u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];
    
    p = [0, 1;
         0, 1;
         0, 2;
         0, 2];
    s = [1, 1, 0, 0, 1, 1, 0, 0;
         1, 1, 0, 0, 1, 1, 0, 0;
         0, 0, 1, 1, 0, 0, 1, 1;
         0, 0, 1, 1, 0, 0, 1, 1];
    u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];
    
    p = [0, 1;
         0, 2;
         1, 1;     
         0, 3];
    s = [1, 0, 0, 0, 1, 0, 0, 0;
         0, 1, 0, 0, 0, 1, 0, 0;
         0, 0, 1, 0, 0, 0, 1, 0;
         0, 0, 0, 1, 0, 0, 0, 1];
    u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

    % Random
%     p = [0, 1;
%          0, 2;
%          0, 3;
%          2, 1];
%     s = [1, 0, 0, 1, 1, 0, 0, 0;
%          0, 1, 0, 0, 0, 1, 0, 0;
%          0, 0, 1, 0, 0, 0, 1, 0;
%          0, 0, 0, 0, 0, 0, 0, 1];
%     u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];
    
    p = [0, 1;
         0, 2;
         0, 3;
         0, 4];
    s = [1, 0, 0, 0, 0, 0, 1, 0;
         0, 0, 0, 1, 0, 1, 0, 0;
         0, 1, 0, 0, 0, 0, 0, 1;
         0, 0, 1, 0, 1, 0, 0, 0];
    u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];
    
%     % Pareto
%     p = [4, 1;
%          3, 1;
%          0, 1;
%          0, 2];
%     s = [1, 0, 0, 0, 1, 0, 0, 0;
%          0, 1, 0, 0, 0, 1, 1, 0;
%          0, 0, 1, 0, 0, 0, 0, 0;
%          0, 0, 0, 1, 0, 0, 0, 1];
%     u_x = [u_x; reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];
end
sys.g = 9.81; % m/s^2
sys.dt = 0.001;
sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(1) = pi;
sys.goal = sys.l_point;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.gamma_ = 0.997;

Op.max_iter = 3000;
Op.max_policy_iter = 300;
Op.gtol = 1e-5;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
Op.save_dir = 'data';
Op.reuse_policy = true;

policies = cell(size(u_x,1), 1);
value = cell(size(u_x,1), 1);
info = cell(size(u_x,1), 1);

for dd=1:1:size(u_x, 1)
    p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    
    disp(sprintf('Decomposition %d/%d', dd, size(u_x,1)));
    sys.decomposition_id = dd;
    [policies{dd,1}, value{dd,1}, info{dd,1}] = dp_decomposition_gpu(sys, Op, p, s);
    policies{dd,1} = cellfun(@(x) gather(x), policies{dd,1}, 'UniformOutput', false);
    value{dd,1} = gather(value{dd,1});
    info{dd,1}.state_grid = cellfun(@(x) gather(x), info{dd,1}.state_grid, 'UniformOutput', false);
end

% p_joint = [zeros(n,1), ones(n,1)];
% s_joint = ones(sys.U_DIMS, sys.X_DIMS);
% sys.decomposition_id = 0;

% disp('Joint');
% [policies_joint, value_joint, info_joint] = dp_decomposition_gpu(sys, Op, p_joint, s_joint);

state_bounds = [repmat([-pi/3, pi/3], [n,1]);
               repmat([-0.5, 0.5], [n,1])];
state_bounds(1,:) = state_bounds(1,:) + pi;
% state_bounds = mat2cell(state_bounds, ones(2*n,1), 2);

% valid_range = cellfun(@(x,y) (x>y(1)) & (x<y(2)), info_joint.state_grid, state_bounds, 'UniformOutput', false);
% valid_range = permute(reshape(permute(cell2mat(valid_range), [linspace(2, 2*n, 2*n-1), 1]), [Op.num_points, 2*n]), [2*n, linspace(1,2*n-1,2*n-1), 2*n+1]);
% valid_range = logical(prod(valid_range, 2*n+1));

% err_dp = zeros(1, size(u_x,1));
% for dd=1:1:size(u_x,1)
%     err_dp(dd) = mean(abs(value{dd,1}(valid_range) - value_joint(valid_range)), 'all');
% end

% save(strcat(Op.save_dir, '/', system_name, '/summary.mat'), 'u_x', 'policies', 'value', 'info', 'policies_joint', 'value_joint', 'info_joint', 'err_dp', 'sys', 'Op');
% save(strcat(Op.save_dir, '/', system_name, '/summary_GA_MCTS_Random.mat'), 'u_x', 'policies', 'value', 'info', 'sys', 'Op', '-v7.3', '-nocompression');

% Closed-loop trajectories
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
% policy_interpolant{size(u_x,1)+1, 1} = cellfun(@(x) griddedInterpolant(info_joint.state_grid{:}, x), ...
%                                                policies_joint{1,3}, 'UniformOutput', false);
% policy_function{size(u_x,1)+1, 1} = cellfun(@(x) @(y) x(y(policies_joint{1,2})), policy_interpolant{size(u_x,1)+1, 1}, 'UniformOutput', false);

% state_bounds = cell2mat(state_bounds);
time_span = [0, 12];
NUM_CTRL = round(time_span(2) / sys.dt);
num_starts = 50;
starts = state_bounds(:,1) + repmat(state_bounds(:,2) - state_bounds(:,1), [1, num_starts]).*rand(sys.X_DIMS, num_starts);
trajectories = cell(size(policy_function, 1), num_starts, 2);
costs = zeros(num_starts, size(policy_function, 1));
sys.X_DIMS_FREE = 1:1:sys.X_DIMS;
for dd=3:1:3
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