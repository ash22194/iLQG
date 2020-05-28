clear;
clc;
close all;
addpath('biped2d');

%% Common parameters
sys.full_DDP = false;   
sys.m = 72;
sys.I = 3;
sys.l0 = 1.05;
sys.g = 9.81;
sys.d = 0.2;
sys.df = 0.5;
sys.T = 5;
sys.dt = 0.001;
NUM_CTRL = round(sys.T / sys.dt);
sys.gamma_ = 0.997;
sys.Q = diag([100, 200, 2, 2, 1000, 10]);
sys.R = 0.0001*eye(4);

lg = 0.96;
alpha1g = pi/2 + asin(sys.df/2/lg);
l2g = sqrt((sys.df + lg*cos(alpha1g))^2 + (lg*sin(alpha1g))^2);
alpha2g = acos((sys.df + lg*cos(alpha1g))/l2g);
sys.goal = [lg; alpha1g; 0; 0; 0; 0];
sys.l_point = sys.goal;
sys.lims = [0, 1.5*sys.m*sys.g;
            0, 1.5*sys.m*sys.g;
            -0.125*sys.m*sys.g, 0.125*sys.m*sys.g;
            -0.125*sys.m*sys.g, 0.125*sys.m*sys.g];
sys.u0 = [sys.m*sys.g*cos(alpha2g)/sin(alpha1g - alpha2g);
          -sys.m*sys.g*cos(alpha1g)/sin(alpha1g - alpha2g);
          0;
          0];

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
Op.Alpha = [0.1];

% Define starts
% com_pos = [0.85, 0.85, 0.9, 0.9, 0.95, 0.95;
%             0.1,  0.4, 0.1, 0.4,  0.1,  0.4];
com_pos = [0.92, 0.92, 1.0, 1.0;
           0.1,  0.4, 0.1, 0.4];
com_pos(2,:) = pi/2 + com_pos(2,:);
com_vel = [ 0.1, -0.1, 0.1, -0.1;
           -0.3, -0.3, -0.4, -0.4];
theta_starts = [-0.3,  -0.15, 0.15, 0.3;
                   0,      0,    0,   0];

x_starts   = nan(6, size(com_pos,2)*size(com_vel,2)*size(theta_starts,2));
com_starts = nan(4, size(com_pos,2)*size(com_vel,2));
ind = 1;
for ii=1:1:size(com_pos, 2)
    for jj=1:1:size(com_vel, 2)
        com_starts(:,(ii-1)*size(com_vel, 2) + jj) = [com_pos(:, ii); com_vel(:, jj)];
        for kk=1:1:size(theta_starts, 2)
            start = [com_pos(:, ii); com_vel(:, jj); theta_starts(:, kk)];
            x_starts(:, ind) = start;
            ind = ind + 1;
        end
    end
end

% alphainit = pi/2 + 0.3;
% x_starts = [sys.df/2/sin(alphainit - pi/2); alphainit; 0; 0; 0; 0];
x_starts = sys.goal;

save_dir = "data/";
save_file = "iLQGBiped2D";

%% Joint
disp('**** Joint ****');
sys_joint = sys;
sys_joint.U_DIMS_FREE = [1;2;3;4];
sys_joint.U_DIMS_FIXED = [];
sys_joint.X_DIMS_FREE = [1;2;3;4;5;6];
sys_joint.X_DIMS_FIXED = [];

Ops.lims = sys_joint.lims;
XJoint = zeros(4, NUM_CTRL+1, size(x_starts, 2));
UJoint = zeros(length(sys_joint.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KJoint = zeros(length(sys_joint.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
TraceJoint = cell(size(x_starts, 2), 1);
timeJoint = zeros(size(x_starts, 2), 1);

for jj=1:1:size(x_starts, 2)
    dyn_joint = @(x, u, i) biped2d_dyn_first_cst(sys_joint, x, u, sys_joint.full_DDP);
    % Uinit = zeros(length(sys_joint.U_DIMS_FREE), NUM_CTRL);
    Uinit = repmat(sys.u0, [1, NUM_CTRL]);
    tic;
    [XJoint(sys_joint.X_DIMS_FREE,:, jj), ...
     UJoint(:,1:NUM_CTRL, jj), ...
     KJoint(:,sys_joint.X_DIMS_FREE,1:NUM_CTRL, jj), TraceJoint{jj,1}] = iLQG(dyn_joint, ...
                                                    x_starts(sys_joint.X_DIMS_FREE, jj), ...
                                                    Uinit, Op);
    timeJoint(jj) = toc;
    XJoint(sys_joint.X_DIMS_FIXED,:, jj) = repmat(sys_joint.l_point(sys_joint.X_DIMS_FIXED), [1, size(XJoint,2)]);
    jj
    XJoint(:,end,jj)
end