clear;
close all;
clc;

%% 

restoredefaultpath;
n = 2;
system_name = sprintf('manipulator%ddof', n);
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/',system_name,'System.mat'));

assert(isfield(sys, 'name') && strcmp(sys.name, system_name), 'Check loaded system!');

sys.X_DIMS = 2*sys.n; % [thi, ... dthi, ...]
sys.U_DIMS = sys.n;   % [taui]
if (n==2)
    sys.m = [2.5; 0.5]/5; % kg
    sys.l = [0.5; 0.25]/2; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8, 8, 0.5, 0.5])/5;
    sys.R = diag(0.003*(Izz(1)./Izz).^2);
    sys.lims = 5*[-Izz/Izz(1), Izz/Izz(1)]; % action limits
    
    Op.num_points = 31 * ones(1, sys.X_DIMS);
    Op.num_action_samples = 15 * ones(1, sys.U_DIMS);
    
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
    sys.lims = [-16, 16; -7.5, 7.5; -1, 1]; % action limits
    
    Op.num_points = 31 * ones(1, sys.X_DIMS);
    Op.num_action_samples = 15 * ones(1, sys.U_DIMS);
    
    % Define decompositions to test
    p = [linspace(0,n-1,n)', ones(n,1)];
    s = repmat(eye(n), [1, 2]);
    u_x = [reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
    
elseif (n==4)
    sys.m = [5.4; 1.8; 0.6; 0.2]; % kg
    sys.l = [0.2; 0.5; 0.25; 0.125]; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8*ones(1,4), 0.2*ones(1,4)])/2;
    sys.R = diag([0.002; 0.004*(Izz(2)./Izz(2:end))]);
    sys.lims = [-24, 24; -15, 15; -7.5, 7.5; -1, 1]; % action limits
    
    Op.num_points = 31 * ones(1, sys.X_DIMS);
    Op.num_action_samples = 15 * ones(1, sys.U_DIMS);
    
    % Define decompositions to test
    p = [linspace(0,n-1,n)', ones(n,1)];
    s = repmat(eye(n), [1, 2]);
    u_x = [reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
end
sys.g = 9.81; % m/s^2
sys.dt = 0.001;
sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(1) = pi;
sys.goal = sys.l_point;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.gamma_ = 0.997;
sys.limits = [0, 2*pi; repmat([-pi, pi], n-1, 1); repmat([-3, 3], n, 1)];

Op.max_iter = 2000;
Op.max_policy_iter = 100;
Op.gtol = 1e-5;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
Op.save_dir = 'data/ddp';

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