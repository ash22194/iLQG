clear;
close all;
clc;

%%

restoredefaultpath;
system_name = 'biped2d';
addpath('general');
addpath(strcat('new_systems/', system_name));
load(strcat('new_systems/', system_name, '/sys.mat'), 'sys');

sys.m = 72;
sys.I = 3;
sys.l0 = 1.05;
sys.g = 9.81;
sys.d = 0.2;
sys.df = 0.5;
sys.T = 4;
sys.dt = 0.001;
sys.X_DIMS = 6;
sys.U_DIMS = 4;

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
sys.fxu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
sys.lims = 2*[0, 1.5*sys.m*sys.g;
              0, 1.5*sys.m*sys.g;
              -0.125*sys.m*sys.g, 0.125*sys.m*sys.g;
              -0.125*sys.m*sys.g, 0.125*sys.m*sys.g];

sys.gamma_ = 0.999;
sys.lambda_ = (1 - sys.gamma_)/sys.dt;
sys.Q = diag([100, 200, 2, 2, 1000, 10]);
sys.R = 0.000002*eye(4);
sys.full_DDP = false;

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
% Op.Alpha = [1];

com_pos = [0.95, 0.95, 1.0, 1.0;
           0.4,  0.3, 0.4, 0.3];
com_pos(2,:) = pi/2 + com_pos(2,:);
com_vel = [ 0.1, -0.1, 0.1, -0.1;
           -0.3, -0.3, 0.3, 0.3];
theta_starts = [-0.2,  -0.2, 0.2,  0.2;
                -0.2,   0.2, -0.2, 0.2];

starts = zeros(sys.X_DIMS, size(com_pos, 2) ...
                           *size(com_vel, 2) ...
                           *size(theta_starts, 2));
count = 0;
for cp=1:1:size(com_pos, 2)
    for cv=1:1:size(com_vel, 2)
        for ts=1:1:size(theta_starts, 2)
            count = count + 1;
            starts(1:2, count) = com_pos(:, cp);
            starts(3:4, count) = com_vel(:, cv);
            starts(5:6, count) = theta_starts(:, ts);
        end
    end
end
starts = starts(:,[1, end]);

%% Test Decompositions
p = [3,1;
     3,1;
     0,1;
     0,1];
s = [1,1,1,1,1,1;
     1,1,1,1,1,1;
     0,0,0,0,0,0;
     0,0,0,0,0,0];
u_x = [reshape(p, 1,2*sys.U_DIMS), reshape(s, 1,sys.U_DIMS*sys.X_DIMS)];

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
sys.name = 'biped2dwu0';
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
sys.name = 'biped2d';
p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjoint, Ujoint, cjoint] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);

sys.u0init = true;
sys.name = 'biped2dwu0';
p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjointwu0, Ujointwu0, cjointwu0] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);