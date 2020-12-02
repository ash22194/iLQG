clear;
close all;
clc;

%% 

restoredefaultpath;
n = 3;
system_name = sprintf('manipulator%ddof', n);
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/',system_name,'System.mat'));

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
    load('data/manipulator3dof/manipulator3dof_paretofront.mat');
    u_x = u_xp;
    
elseif (n==4)
    sys.m = [5.4; 1.8; 0.6; 0.2]; % kg
    sys.l = [0.2; 0.5; 0.25; 0.125]; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8*ones(1,4), 0.2*ones(1,4)])/2;
    sys.R = diag([0.002; 0.004*(Izz(2)./Izz(2:end))]);
    sys.limits = [pi/2, 3*pi/2; repmat([-pi/2, pi/2], n-1, 1); repmat([-1.5, 1.5], n, 1)];
    sys.lims = [-24, 24; -15, 15; -7.5, 7.5; -1, 1]; % action limits
    
    Op.num_points = [9, 9, 9, 9, 7, 7, 7, 7];
    Op.num_action_samples = [12, 7, 3, 2];
    
    % Define decompositions to test
    load('data/manipulator4dof/manipulator4dof_paretofront.mat');
    u_x = u_xp;
end
sys.g = 9.81; % m/s^2
sys.dt = 0.001;
sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(1) = pi;
sys.goal = sys.l_point;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.gamma_ = 0.997;

Op.max_iter = 2000;
Op.max_policy_iter = 100;
Op.gtol = 1e-5;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
Op.save_dir = 'data';
Op.reuse_policy = true;

policies = cell(size(u_x,1), 1);
value = cell(size(u_x,1), 1);
info = cell(size(u_x,1), 1);

for dd=1:1:size(u_x,1)
    p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    
    disp(sprintf('Decomposition %d/%d', dd, size(u_x,1)));
    sys.decomposition_id = dd;
    [policies{dd,1}, value{dd,1}, info{dd,1}] = dp_decomposition(sys, Op, p, s);
end

p_joint = [zeros(n,1), ones(n,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
disp('Joint');
sys.decomposition_id = 0;
[policies_joint, value_joint, info_joint] = dp_decomposition(sys, Op, p_joint, s_joint);

state_bounds = [repmat([-pi/3, pi/3], [n,1]);
               repmat([-0.5, 0.5], [n,1])];
state_bounds(1,:) = state_bounds(1,:) + pi;
state_bounds = mat2cell(state_bounds, ones(2*n,1), 2);

valid_range = cellfun(@(x,y) (x>y(1)) & (x<y(2)), info_joint.state_grid, state_bounds, 'UniformOutput', false);
valid_range = permute(reshape(permute(cell2mat(valid_range), [linspace(2, 2*n, 2*n-1), 1]), [Op.num_points, 2*n]), [2*n, linspace(1,2*n-1,2*n-1), 2*n+1]);
valid_range = logical(prod(valid_range, 2*n+1));

err_dp = zeros(1, size(u_x,1));
for dd=1:1:size(u_x,1)
    err_dp(dd) = mean(abs(value{dd,1}(valid_range) - value_joint(valid_range)), 'all');
end

save(strcat(Op.save_dir, '/', system_name, '/summary.mat'), 'u_x', 'policies', 'value', 'info', 'policies_joint', 'value_joint', 'info_joint', 'err_dp', 'sys', 'Op');
