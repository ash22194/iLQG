clear;
close all;
clc;

%% 

restoredefaultpath;
system_name = 'manipulator3dof';
policy_type = 'cascaded';
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/',system_name,'System.mat'));

assert(isfield(sys, 'name') && strcmp(sys.name, system_name), 'Check loaded system!');

if (strcmp(system_name, 'cartpole'))
    sys.mc = 5;
    sys.mp = 1;
    sys.l = 0.9;
    sys.g = 9.81;
    sys.dt = 0.001;
    sys.X_DIMS = 4; % x, x_dot, th, th_dot
    sys.U_DIMS = 2; % F, tau 
    
    sys.goal = [0;0;pi;0];
    sys.u0 = zeros(sys.U_DIMS,1);
    sys.Q = diag([25, 0.02, 25, 0.02]);
    sys.R = 0.001*eye(sys.U_DIMS);
    sys.gamma_ = 0.997;
    sys.limits = [-1.5, 1.5; -3, 3; 0, 2*pi; -3, 3];
    sys.lims   = 9*[-ones(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
    
    Op.num_points = 31 * ones(1,sys.X_DIMS);
    Op.num_action_samples = 15 * ones(1, sys.U_DIMS);
    Op.max_iter = 2000;
    Op.max_policy_iter = 100;
    Op.gtol = 1e-5;
    Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
    Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
    
    if (strcmp(policy_type, 'joint'))
        p = [0, 1;
             0, 1];
        s = [1,1,1,1;
             1,1,1,1];
    elseif (strcmp(policy_type, 'cascaded'))
        p = [0, 1;
             1, 1];
        s = [1,1,0,0;
             0,0,1,1];
    elseif (strcmp(policy_type, 'decoupled'))
        p = [0, 1;
             0, 2];
        s = [1,1,0,0;
             0,0,1,1];
    end

    
elseif (strcmp(system_name, 'biped2d'))
    sys.m = 72;
    sys.I = 3;
    sys.l0 = 1.05;
    sys.d = 0.2;
    sys.df = 0.5;
    sys.g = 9.81;
    sys.dt = 0.001;
    sys.X_DIMS = 6; % l, alpha, x_dot, z_dot, th, th_dot
    sys.U_DIMS = 4; % F1, F2, tau1, tau2
    
    lg = 0.96;
    alpha1g = pi/2 + asin(sys.df/2/lg);
    l2g = sqrt((sys.df + lg*cos(alpha1g))^2 + (lg*sin(alpha1g))^2);
    alpha2g = acos((sys.df + lg*cos(alpha1g))/l2g);
    sys.goal = [lg; alpha1g; 0; 0; 0; 0];
    sys.u0 = [m*g*cos(alpha2g)/sin(alpha1g - alpha2g); -m*g*cos(alpha1g)/sin(alpha1g - alpha2g); 0; 0];
    sys.Q = diag([100, 200, 2, 2, 1000, 10]);
    sys.R = 0.000002*eye(4);
    sys.gamma_ = 0.999;
    sys.limits = [sys.l0 - 0.6, sys.l0 + 0.1;
                  pi/2, pi/2 + 0.6;
                  -0.3, 0.3;
                  -0.4, 0.4;
                  -pi/6, pi/6;
                  -1.5, 1.5];
    sys.lims = 2*[0, 1.5*sys.m*sys.g;
                  0, 1.5*sys.m*sys.g;
                  -0.125*sys.m*sys.g, 0.125*sys.m*sys.g;
                  -0.125*sys.m*sys.g, 0.125*sys.m*sys.g];
    
    Op.num_points = [14, 16, 17, 19, 21, 21];
    Op.num_action_samples = [15, 15, 5, 5];
    Op.max_iter = 1000;
    Op.max_policy_iter = 100;
    Op.gtol = 2e-6;
    Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
    Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
    
    if (strcmp(policy_type, 'joint'))
        p = [zeros(4,1), ones(4,1)];
        s = ones(sys.U_DIMS, sys.X_DIMS);
    elseif (strcmp(policy_type, 'cascaded'))
        p = [3, 1; 3, 1; 0, 1; 0, 1];
        s = [1,1,1,1,0,0;
             1,1,1,1,0,0;
             0,0,0,0,1,1;
             0,0,0,0,1,1];
    elseif (strcmp(policy_type, 'decoupled'))
        p = [0, 1; 0, 1; 0, 2; 0, 2];
        s = [1,1,1,1,0,0;
             1,1,1,1,0,0;
             0,0,0,0,1,1;
             0,0,0,0,1,1];
    end

elseif (strcmp(system_name, 'manipulator2dof'))
    sys.m = [2.5; 0.5]/5; % kg
    sys.l = [0.5; 0.25]/2; % m
    Izz = sys.m.*((sys.l));
    sys.g = 9.81;
    sys.dt = 0.001;
    sys.X_DIMS = 4; % th1, th2, dth1, dth2
    sys.U_DIMS = 2; % tau1, tau2
    
    sys.goal = [pi;0;0;0];
    sys.u0 = zeros(sys.U_DIMS,1);
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8, 8, 0.5, 0.5]);
    sys.R = diag(0.003*(Izz(1)./Izz).^2);
    sys.gamma_ = 0.997;
    sys.limits = [0, 2*pi; -pi, pi; -3, 3; -3, 3];
    sys.lims = 5*[-Izz/Izz(1), Izz/Izz(1)]; % action limits
    
    Op.num_points = 31 * ones(1, sys.X_DIMS);
    Op.num_action_samples = [5, 5];
    Op.max_iter = 2000;
    Op.max_policy_iter = 100;
    Op.gtol = 1e-5;
    Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
    Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
    
    if (strcmp(policy_type, 'joint'))
        p = [0, 1;
             0, 1];
        s = [1,1,1,1;
             1,1,1,1];
    elseif (strcmp(policy_type, 'cascaded'))
        p = [0, 1;
             1, 1];
        s = [1,0,1,0;
             0,1,0,1];
    elseif (strcmp(policy_type, 'decoupled'))
        p = [0, 1;
             0, 2];
        s = [1,0,1,0;
             0,1,0,1];
    end
    
elseif (strcmp(system_name, 'manipulator3dof'))
    sys.m = [2.75; 0.55; 0.11]; % kg
    sys.l = [0.5; 0.25; 0.125]; % m
    sys.g = 9.81;
    sys.dt = 0.001;
    sys.X_DIMS = 6; % th1, th2, th3, dth1, dth2, dth3
    sys.U_DIMS = 3; % tau1, tau2, tau3
    
    sys.goal = [pi;0;0;0;0;0];
    sys.u0 = zeros(sys.U_DIMS,1);
    Izz = sys.m.*((sys.l));
    sys.Q = diag([1.6*ones(1,3), 0.12*ones(1,3)]);
    sys.R = diag(0.004*(Izz(1)./Izz));
    sys.gamma_ = 0.997;
    sys.limits = [0, 2*pi; 0, 2*pi; 0, 2*pi; -3, 3; -3, 3; -3, 3];
    sys.lims = [-16, 16; -7.5, 7.5; -1, 1]; % action limits
    
    Op.num_points = 21 * ones(1, sys.X_DIMS);
    Op.num_action_samples = [15, 10, 5];
    Op.max_iter = 2000;
    Op.max_policy_iter = 100;
    Op.gtol = 1e-5;
    Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
    Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
    
    if (strcmp(policy_type, 'joint'))
        p = [zeros(3,1), ones(3,1)];
        s = ones(sys.U_DIMS, sys.X_DIMS);
    elseif (strcmp(policy_type, 'cascaded'))
        p = [linspace(0,2,3)', ones(3,1)];
        s = repmat(eye(3), [1,2]);
    elseif (strcmp(policy_type, 'decoupled'))
        p = [zeros(3,1), linspace(1,3,3)'];
        s = repmat(eye(3), [1,2]);
    end
    
end

% profile on;
policies = dp_decomposition(sys, Op, p, s);
% profile off;
% profsave(profile('info'), 'data/DPcodeprofile');

save(strcat('data/', system_name, '_', policy_type, '_policies.mat'), 'policies', 'sys');