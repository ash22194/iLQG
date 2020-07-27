clear;
clc;
close all;
addpath('biped2d');
addpath('iLQG utilities/kdtree/kdtree/lib');

%% Common parameters
sys.full_DDP = false;   
sys.m = 72;
sys.I = 3;
sys.l0 = 1.05;
sys.g = 9.81;
sys.d = 0.2;
sys.df = 0.5;
sys.T = 4;
sys.dt = 0.001;
NUM_CTRL = round(sys.T / sys.dt);
sys.gamma_ = 0.999;
sys.Q = diag([100, 200, 2, 2, 1000, 10]);
sys.R = 0.000002*eye(4);

lg = 0.96;
alpha1g = pi/2 + asin(sys.df/2/lg);
l2g = sqrt((sys.df + lg*cos(alpha1g))^2 + (lg*sin(alpha1g))^2);
alpha2g = acos((sys.df + lg*cos(alpha1g))/l2g);
sys.goal = [lg; alpha1g; 0; 0; 0; 0];
sys.l_point = sys.goal;
sys.lims = 2*[0, 1.5*sys.m*sys.g;
              0, 1.5*sys.m*sys.g;
              -0.125*sys.m*sys.g, 0.125*sys.m*sys.g;
              -0.125*sys.m*sys.g, 0.125*sys.m*sys.g];
sys.u0 = [sys.m*sys.g*cos(alpha2g)/sin(alpha1g - alpha2g);
          -sys.m*sys.g*cos(alpha1g)/sin(alpha1g - alpha2g);
          0;
          0];

%% Linearized system
x = sym('x', [6, 1]); % l, alpha1, x_dot, z_dot, th, th_dot
u = sym('u', [4, 1]); % F_leg1, F_leg2, T_hip1, T_hip2
x_hip     = x(1)*cos(x(2));
z_hip     = x(1)*sin(x(2));
l2        = sqrt((x_hip + sys.df)^2 + z_hip^2);
a2        = acos((sys.df + x_hip)/l2);
x_hip_dot = x(3) + sys.d*cos(x(5))*x(6);
z_hip_dot = x(4) + sys.d*sin(x(5))*x(6);

zero_dyn = [0; 0; 0; -sys.g; 0; 0];
state_dyn = [x_hip_dot*cos(x(2)) + z_hip_dot*sin(x(2));
             (-x_hip_dot*sin(x(2)) + z_hip_dot*cos(x(2)))/x(1);
             0;
             0;
             x(6);
             0];
act_dyn = [0;
           0;
           (u(1)*cos(x(2)) + u(2)*cos(a2) ...
            + u(3)/x(1)*sin(x(2)) + u(4)/l2*sin(a2))/sys.m;
           (u(1)*sin(x(2)) + u(2)*sin(a2) ...
            - u(3)/x(1)*cos(x(2)) - u(4)/l2*cos(a2))/sys.m;
           0;
           (u(3)*(1 + sys.d/x(1)*sin(x(2) - x(6))) + u(4)*(1 + sys.d/l2*sin(a2 - x(6))) ...
            + u(1)*sys.d*cos(x(2) - x(6)) + u(2)*sys.d*cos(a2 - x(6)))/sys.I];
fxu = zero_dyn + state_dyn + act_dyn;
state_dyn_x = jacobian(fxu,x);
act_dyn_u = jacobian(fxu,u);
A = eval(subs(state_dyn_x, [x; u], [sys.l_point; sys.u0])); 
B = eval(subs(act_dyn_u, [x; u], [sys.l_point; sys.u0]));
lambda_ = (1 - sys.gamma_)/sys.dt;

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
% Op.Alpha = [1];

% Define starts
com_pos = [0.95, 0.95, 1.0, 1.0;
           0.4,  0.3, 0.4, 0.3];
com_pos(2,:) = pi/2 + com_pos(2,:);
com_vel = [ 0.1, -0.1, 0.1, -0.1;
           -0.3, -0.3, 0.3, 0.3];
theta_starts = [-0.2,  -0.2, 0.2,  0.2;
                -0.2,   0.2, -0.2, 0.2];

x_starts = nan(6, size(com_pos,2)*size(com_vel,2)*size(theta_starts,2));
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

% x_starts = x_starts(:,1);
% com_starts = com_starts(:,1);
% theta_starts = theta_starts(:,1);

save_dir = "data/";
save_file = strcat("iLQGBiped2DDecomposed_cornerstart_all8");

%% Joint
[K_joint, S_joint, e_joint] = lqr(A - lambda_/2*eye(size(A,1)), B, ...
                                  sys.Q, sys.R, ...
                                  zeros(size(A,1), size(B,2)));

sys_joint = sys;
sys_joint.U_DIMS_FREE = [1;2;3;4];
sys_joint.U_DIMS_FIXED = [];
sys_joint.X_DIMS_FREE = [1;2;3;4;5;6];
sys_joint.X_DIMS_FIXED = [];

Ops.lims = sys_joint.lims;
XJoint = zeros(4, NUM_CTRL+1, size(x_starts, 2));
UJoint = zeros(length(sys_joint.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KJoint = zeros(length(sys_joint.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
XJointinit = nan(6, NUM_CTRL+1, size(x_starts, 2));
UJointinit = nan(4, NUM_CTRL+1, size(x_starts, 2));
CostJointinit = zeros(size(x_starts, 2), 1);
TraceJoint = cell(size(x_starts, 2), 1);
timeJoint = zeros(size(x_starts, 2), 1);

for jj=1:1:size(x_starts, 2)
    dyn_joint = @(x, u, i) biped2d_dyn_first_cst(sys_joint, x, u, sys_joint.full_DDP);
    
    XJointinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UJointinit(:,ii, jj) = min(max(sys_joint.u0 - K_joint*(XJointinit(:,ii, jj) - sys_joint.goal), sys_joint.lims(:,1)), sys_joint.lims(:,2));
        XJointinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_joint, XJointinit(:,ii, jj), UJointinit(:,ii, jj), sys_joint.dt);
        if (any(isnan(UJointinit(:,ii, jj))) || any(isnan(XJointinit(:,ii, jj))))
            disp('Check X, U Joint');
        end
        CostJointinit(jj) = CostJointinit(jj) + discount*l_Biped2DFirst(sys_joint, XJointinit(:,ii, jj), UJointinit(:,ii, jj) - sys_joint.u0)*sys_joint.dt;
        discount = discount*sys_joint.gamma_;
    end
    CostJointinit(jj) = CostJointinit(jj) + discount*l_Biped2DFirst(sys_joint, XJointinit(:,NUM_CTRL+1, jj), zeros(4,1))*sys_joint.dt;
    
    Op.cost = CostJointinit(jj);
    tic;
    [XJoint(sys_joint.X_DIMS_FREE,:, jj), ...
     UJoint(:,1:NUM_CTRL, jj), ...
     KJoint(:,sys_joint.X_DIMS_FREE,1:NUM_CTRL, jj), TraceJoint{jj,1}] = iLQG(dyn_joint, ...
                                                    XJointinit(:,:,jj), ...
                                                    UJointinit(:,1:NUM_CTRL,jj), Op);
    timeJoint(jj) = toc;
    XJoint(sys_joint.X_DIMS_FIXED,:, jj) = repmat(sys_joint.l_point(sys_joint.X_DIMS_FIXED), [1, size(XJoint,2)]);
    jj
end

%% Full - F, Full - T
disp('**** F - Full, T - Full ****');
% F Full
% A_ = A;
% B_ = B(:,1:2);

disp('F - Full')
sys_FFull = sys;
sys_FFull.U_DIMS_FREE = [1; 2];
sys_FFull.U_DIMS_FIXED = linspace(1,4,4)';
sys_FFull.U_DIMS_FIXED(sys_FFull.U_DIMS_FREE) = [];
sys_FFull.X_DIMS_FREE = [1;2;3;4;5;6];
sys_FFull.X_DIMS_FIXED = [];

u0 = sys.u0;
u0(sys_FFull.U_DIMS_FIXED) = 0;
A_ = eval(subs(state_dyn_x, [x; u], [sys.l_point; u0]));
B_ = eval(subs(act_dyn_u, [x; u], [sys.l_point; u0]));
B_ = B_(:,1:2);
Q_ = sys.Q;
R_ = sys.R(1:2,1:2);
[K_FFull, S_FFull, e_FFull] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);

Op.lims = sys_FFull.lims(sys_FFull.U_DIMS_FREE, :);
XFFull = zeros(length(sys_FFull.X_DIMS_FREE) + length(sys_FFull.X_DIMS_FIXED), NUM_CTRL+1, size(x_starts, 2));
UFFull = zeros(length(sys_FFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KFFull = zeros(length(sys_FFull.U_DIMS_FREE), length(sys_FFull.X_DIMS_FREE) + length(sys_FFull.X_DIMS_FIXED), ...
               NUM_CTRL+1, size(x_starts, 2));
XFFullinit = nan(length(sys_FFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UFFullinit = nan(length(sys_FFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostFFullinit = zeros(size(x_starts, 2), 1);
TraceFFull = cell(size(x_starts, 2), 1);
timeFFull = zeros(size(x_starts, 2), 1);

u0 = sys_FFull.u0(sys_FFull.U_DIMS_FREE, 1);
goal = sys_FFull.goal(sys_FFull.X_DIMS_FREE, 1);
lims = sys_FFull.lims(sys_FFull.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_FFull = @(x, u, i) biped2d_dyn_first_cst(sys_FFull, x, u, sys_FFull.full_DDP);
    
    XFFullinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UFFullinit(:,ii, jj) = min(max(u0 - K_FFull*(XFFullinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XFFullinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_FFull, XFFullinit(:,ii, jj), UFFullinit(:,ii, jj), sys_FFull.dt);
        if (any(isnan(UFFullinit(:,ii, jj))) || any(isnan(XFFullinit(:,ii, jj))))
            disp('Check X, U FFull');
        end
        CostFFullinit(jj) = CostFFullinit(jj) + discount*l_Biped2DFirst(sys_FFull, XFFullinit(:,ii, jj), UFFullinit(:,ii, jj) - u0)*sys_FFull.dt;
        discount = discount*sys_FFull.gamma_;
    end
    CostFFullinit(jj) = CostFFullinit(jj) + discount*l_Biped2DFirst(sys_FFull, XFFullinit(:,NUM_CTRL+1, jj), zeros(2,1))*sys_FFull.dt;
    
    Op.cost = CostFFullinit(jj);
    tic;
    [XFFull(sys_FFull.X_DIMS_FREE,:, jj), ...
     UFFull(:,1:NUM_CTRL, jj), ...
     KFFull(:,sys_FFull.X_DIMS_FREE,1:NUM_CTRL, jj), TraceFFull{jj,1}] = iLQG(dyn_FFull, ...
                                                    XFFullinit(:,:, jj), ...
                                                    UFFullinit(:,1:NUM_CTRL, jj), Op);
    timeFFull(jj) = toc;
    jj
end

% T Full
A_ = [A - B(:,1:2)*K_FFull];
B_ = [zeros(6,2), B(:,3:4)];
Q_ = sys.Q + K_FFull'*sys.R(1:2,1:2)*K_FFull;
R_ = sys.R;
[K_TSFull, S_TSFull, e_TSFull] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);
K_TSFull = K_TSFull(3:4,:);

disp('T - Both')
sys_TSFull = sys;
sys_TSFull.U_DIMS_FREE = [3; 4];
sys_TSFull.U_DIMS_FIXED = linspace(1,4,4)';
sys_TSFull.U_DIMS_FIXED(sys_TSFull.U_DIMS_FREE) = [];
sys_TSFull.X_DIMS_FREE = [1;2;3;4;5;6];
sys_TSFull.X_DIMS_FIXED = [];

Op.lims = sys_TSFull.lims(sys_TSFull.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [1;2;3;4;5;6];
XTSFull = zeros(length(sys_TSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTSFull = zeros(length(sys_TSFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTSFull = zeros(length(sys_TSFull.U_DIMS_FREE), length(sys_TSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XFFullClose = nan(length(sys_TSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UFFullClose = zeros(length(sys_FFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KFFullClose = zeros(length(sys_FFull.U_DIMS_FREE), length(sys_TSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTSFullinit = nan(length(sys_TSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTSFullinit = nan(length(sys_TSFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostTSFullinit = zeros(size(x_starts, 2), 1);
TraceTSFull = cell(size(x_starts, 2), 1);
timeTSFull = zeros(size(x_starts, 2), 1);

u0Free = sys_TSFull.u0(sys_TSFull.U_DIMS_FREE, 1);
u0Fixed = sys_TSFull.u0(sys_TSFull.U_DIMS_FIXED, 1);
goal = sys_TSFull.goal(sys_TSFull.X_DIMS_FREE, 1);
lims = sys_TSFull.lims(sys_TSFull.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_TSFull = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_TSFull, x, u, k, K, xn, sys_TSFull.full_DDP);
    
    XTSFullinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UTSFullinit(:,ii, jj) = min(max(u0Free - K_TSFull*(XTSFullinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XTSFullinit(:,ii+1, jj) = f_Biped2DSecond_finite(sys_TSFull, XTSFullinit(:,ii, jj), UTSFullinit(:,ii, jj), u0Fixed, -K_FFull, sys_FFull.goal, sys_TSFull.dt);
        if (any(isnan(UTSFullinit(:,ii, jj))) || any(isnan(XTSFullinit(:,ii, jj))))
            disp('Check X, U TSFull');
        end
        CostTSFullinit(jj) = CostTSFullinit(jj) + discount*l_Biped2DSecond(sys_TSFull, XTSFullinit(:,ii, jj), UTSFullinit(:,ii, jj) - u0Free, zeros(length(u0Fixed),1), -K_FFull, sys_FFull.goal)*sys_TSFull.dt;
        discount = discount*sys_TSFull.gamma_;
    end
    CostTSFullinit(jj) = CostTSFullinit(jj) + discount*l_Biped2DSecond(sys_TSFull, XTSFullinit(:,NUM_CTRL+1, jj), zeros(length(u0Free),1), zeros(length(u0Fixed),1), -K_FFull, sys_FFull.goal)*sys_TSFull.dt;
    
    Op.cost = CostTSFullinit(jj);
    tic;
    [XTSFull(:,:, jj), ...
     UTSFull(:,1:NUM_CTRL, jj), ...
     KTSFull(:,:,1:NUM_CTRL, jj), ...
     UFFullClose(:,:, jj), ...
     KFFullClose(:,:,:, jj), ...
     XFFullClose(:,:, jj), TraceTSFull{jj,1}] = iLQGSecondKDTree(dyn_TSFull, ...
                                        XTSFullinit(:,:, jj), ...
                                        UTSFullinit(:,1:NUM_CTRL, jj), ...
                                        UFFull(:,:, jj), ...
                                        KFFull(:,:,:, jj), ...
                                        XFFull(sys_FFull.X_DIMS_FREE,:, jj), ...
                                        Op);
    timeTSFull(jj) = toc;
    jj

end

%% Full - T, Full - F
disp('**** T - Full, F - Full ****');
% T Full
% A_ = A;
% B_ = B(:,3:4);

disp('T - Full')
sys_TFull = sys;
sys_TFull.U_DIMS_FREE = [3; 4];
sys_TFull.U_DIMS_FIXED = linspace(1,4,4)';
sys_TFull.U_DIMS_FIXED(sys_TFull.U_DIMS_FREE) = [];
sys_TFull.X_DIMS_FREE = [1;2;3;4;5;6];
sys_TFull.X_DIMS_FIXED = [];

u0 = sys.u0;
u0(sys_TFull.U_DIMS_FIXED) = 0;
A_ = eval(subs(state_dyn_x, [x; u], [sys.l_point; u0]));
B_ = eval(subs(act_dyn_u, [x; u], [sys.l_point; u0]));
B_ = B_(:,3:4);
Q_ = sys.Q;
R_ = sys.R(3:4,3:4);
[K_TFull, S_TFull, e_TFull] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);

Op.lims = sys_TFull.lims(sys_TFull.U_DIMS_FREE, :);
XTFull = zeros(length(sys_TFull.X_DIMS_FREE) + length(sys_TFull.X_DIMS_FIXED), NUM_CTRL+1, size(x_starts, 2));
UTFull = zeros(length(sys_TFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTFull = zeros(length(sys_TFull.U_DIMS_FREE), length(sys_TFull.X_DIMS_FREE) + length(sys_FFull.X_DIMS_FIXED), ...
               NUM_CTRL+1, size(x_starts, 2));
XTFullinit = nan(length(sys_TFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTFullinit = nan(length(sys_TFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostTFullinit = zeros(size(x_starts, 2), 1);
TraceTFull = cell(size(x_starts, 2), 1);
timeTFull = zeros(size(x_starts, 2), 1);

u0 = sys_TFull.u0(sys_TFull.U_DIMS_FREE, 1);
goal = sys_TFull.goal(sys_TFull.X_DIMS_FREE, 1);
lims = sys_TFull.lims(sys_TFull.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_TFull = @(x, u, i) biped2d_dyn_first_cst(sys_TFull, x, u, sys_TFull.full_DDP);
    
    XTFullinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UTFullinit(:,ii, jj) = min(max(u0 - K_TFull*(XTFullinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XTFullinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_TFull, XTFullinit(:,ii, jj), UTFullinit(:,ii, jj), sys_TFull.dt);
        if (any(isnan(UTFullinit(:,ii, jj))) || any(isnan(XTFullinit(:,ii, jj))))
            disp('Check X, U TFull');
        end
        CostTFullinit(jj) = CostTFullinit(jj) + discount*l_Biped2DFirst(sys_TFull, XTFullinit(:,ii, jj), UTFullinit(:,ii, jj) - u0)*sys_TFull.dt;
        discount = discount*sys_TFull.gamma_;
    end
    CostTFullinit(jj) = CostTFullinit(jj) + discount*l_Biped2DFirst(sys_TFull, XTFullinit(:,NUM_CTRL+1, jj), zeros(2,1))*sys_TFull.dt;
    
    Op.cost = CostTFullinit(jj);
    tic;
    [XTFull(sys_TFull.X_DIMS_FREE,:, jj), ...
     UTFull(:,1:NUM_CTRL, jj), ...
     KTFull(:,sys_TFull.X_DIMS_FREE,1:NUM_CTRL, jj), TraceTFull{jj,1}] = iLQG(dyn_TFull, ...
                                                    XTFullinit(:,:, jj), ...
                                                    UTFullinit(:,1:NUM_CTRL, jj), Op);
    timeTFull(jj) = toc;
    jj
end

% F Full
A_ = [A - B(:,3:4)*K_TFull];
B_ = [ B(:,1:2), zeros(6,2)];
Q_ = sys.Q + K_TFull'*sys.R(3:4,3:4)*K_TFull;
R_ = sys.R;
[K_FSFull, S_FSFull, e_FSFull] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);
K_FSFull = K_FSFull(1:2,:);

disp('F - Full')
sys_FSFull = sys;
sys_FSFull.U_DIMS_FREE = [1; 2];
sys_FSFull.U_DIMS_FIXED = linspace(1,4,4)';
sys_FSFull.U_DIMS_FIXED(sys_FSFull.U_DIMS_FREE) = [];
sys_FSFull.X_DIMS_FREE = [1;2;3;4;5;6];
sys_FSFull.X_DIMS_FIXED = [];

Op.lims = sys_FSFull.lims(sys_FSFull.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [1;2;3;4;5;6];
XFSFull = zeros(length(sys_FSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UFSFull = zeros(length(sys_FSFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KFSFull = zeros(length(sys_FSFull.U_DIMS_FREE), length(sys_FSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTFullClose = nan(length(sys_FSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTFullClose = zeros(length(sys_TFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTFullClose = zeros(length(sys_TFull.U_DIMS_FREE), length(sys_FSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XFSFullinit = nan(length(sys_FSFull.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UFSFullinit = nan(length(sys_FSFull.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostFSFullinit = zeros(size(x_starts, 2), 1);
TraceFSFull = cell(size(x_starts, 2), 1);
timeFSFull = zeros(size(x_starts, 2), 1);

u0Free = sys_FSFull.u0(sys_FSFull.U_DIMS_FREE, 1);
u0Fixed = sys_FSFull.u0(sys_FSFull.U_DIMS_FIXED, 1);
goal = sys_FSFull.goal(sys_FSFull.X_DIMS_FREE, 1);
lims = sys_FSFull.lims(sys_FSFull.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_FSFull = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_FSFull, x, u, k, K, xn, sys_FSFull.full_DDP);
    
    XFSFullinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UFSFullinit(:,ii, jj) = min(max(u0Free - K_FSFull*(XFSFullinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XFSFullinit(:,ii+1, jj) = f_Biped2DSecond_finite(sys_FSFull, XFSFullinit(:,ii, jj), UFSFullinit(:,ii, jj), u0Fixed, -K_TFull, sys_TFull.goal, sys_FSFull.dt);
        if (any(isnan(UFSFullinit(:,ii, jj))) || any(isnan(XFSFullinit(:,ii, jj))))
            disp('Check X, U FSFull');
        end
        CostFSFullinit(jj) = CostFSFullinit(jj) + discount*l_Biped2DSecond(sys_FSFull, XFSFullinit(:,ii, jj), UFSFullinit(:,ii, jj) - u0Free, zeros(length(u0Fixed),1), -K_TFull, sys_TFull.goal)*sys_FSFull.dt;
        discount = discount*sys_FSFull.gamma_;
    end
    CostFSFullinit(jj) = CostFSFullinit(jj) + discount*l_Biped2DSecond(sys_FSFull, XFSFullinit(:,NUM_CTRL+1, jj), zeros(length(u0Free),1), zeros(length(u0Fixed),1), -K_TFull, sys_TFull.goal)*sys_FSFull.dt;
    
    Op.cost = CostFSFullinit(jj);
    tic;
    [XFSFull(:,:, jj), ...
     UFSFull(:,1:NUM_CTRL, jj), ...
     KFSFull(:,:,1:NUM_CTRL, jj), ...
     UTFullClose(:,:, jj), ...
     KTFullClose(:,:,:, jj), ...
     XTFullClose(:,:, jj), TraceFSFull{jj,1}] = iLQGSecondKDTree(dyn_FSFull, ...
                                        XFSFullinit(:,:, jj), ...
                                        UFSFullinit(:,1:NUM_CTRL, jj), ...
                                        UTFull(:,:, jj), ...
                                        KTFull(:,:,:, jj), ...
                                        XTFull(sys_TFull.X_DIMS_FREE,:, jj), ...
                                        Op);
    timeFSFull(jj) = toc;
    jj

end


%% COM - F, Torso - T
disp('**** F - COM, T - Both ****');
% COM First
% A_ = A(1:4,1:4);
% B_ = B(1:4,1:2);

disp('F - COM')
sys_COMFF = sys;
sys_COMFF.U_DIMS_FREE = [1; 2];
sys_COMFF.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMFF.U_DIMS_FIXED(sys_COMFF.U_DIMS_FREE) = [];
sys_COMFF.X_DIMS_FREE = [1;2;3;4];
sys_COMFF.X_DIMS_FIXED = linspace(1,6,6)';
sys_COMFF.X_DIMS_FIXED(sys_COMFF.X_DIMS_FREE) = [];

u0 = sys.u0;
u0(sys_COMFF.U_DIMS_FIXED) = 0;
A_ = eval(subs(state_dyn_x, [x; u], [sys.l_point; u0]));
A_ = A_(1:4,1:4);
B_ = eval(subs(act_dyn_u, [x; u], [sys.l_point; u0]));
B_ = B_(1:4,1:2);
Q_ = sys.Q(1:4,1:4);
R_ = sys.R(1:2,1:2);
[K_COMFF, S_COMFF, e_COMFF] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);

Op.lims = sys_COMFF.lims(sys_COMFF.U_DIMS_FREE, :);
XCOMFF = zeros(length(sys_COMFF.X_DIMS_FREE) + length(sys_COMFF.X_DIMS_FIXED), NUM_CTRL+1, size(com_starts, 2));
UCOMFF = zeros(length(sys_COMFF.U_DIMS_FREE), NUM_CTRL+1, size(com_starts, 2));
KCOMFF = zeros(length(sys_COMFF.U_DIMS_FREE), length(sys_COMFF.X_DIMS_FREE) + length(sys_COMFF.X_DIMS_FIXED), ...
               NUM_CTRL+1, size(com_starts, 2));
XCOMFFinit = nan(length(sys_COMFF.X_DIMS_FREE), NUM_CTRL+1, size(com_starts, 2));
UCOMFFinit = nan(length(sys_COMFF.U_DIMS_FREE), NUM_CTRL+1, size(com_starts, 2));
CostCOMFFinit = zeros(size(com_starts, 2), 1);
TraceCOMFF = cell(size(com_starts, 2), 1);
timeCOMFF = zeros(size(com_starts, 2), 1);

u0 = sys_COMFF.u0(sys_COMFF.U_DIMS_FREE, 1);
goal = sys_COMFF.goal(sys_COMFF.X_DIMS_FREE, 1);
lims = sys_COMFF.lims(sys_COMFF.U_DIMS_FREE, :);
for jj=1:1:size(com_starts, 2)
    dyn_COMFF = @(x, u, i) biped2d_dyn_first_cst(sys_COMFF, x, u, sys_COMFF.full_DDP);
    
    XCOMFFinit(:,1, jj) = com_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UCOMFFinit(:,ii, jj) = min(max(u0 - K_COMFF*(XCOMFFinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XCOMFFinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_COMFF, XCOMFFinit(:,ii, jj), UCOMFFinit(:,ii, jj), sys_COMFF.dt);
        if (any(isnan(UCOMFFinit(:,ii, jj))) || any(isnan(XCOMFFinit(:,ii, jj))))
            disp('Check X, U COMFF');
        end
        CostCOMFFinit(jj) = CostCOMFFinit(jj) + discount*l_Biped2DFirst(sys_COMFF, XCOMFFinit(:,ii, jj), UCOMFFinit(:,ii, jj) - u0)*sys_COMFF.dt;
        discount = discount*sys_COMFF.gamma_;
    end
    CostCOMFFinit(jj) = CostCOMFFinit(jj) + discount*l_Biped2DFirst(sys_COMFF, XCOMFFinit(:,NUM_CTRL+1, jj), zeros(2,1))*sys_COMFF.dt;
    
    Op.cost = CostCOMFFinit(jj);
    tic;
    [XCOMFF(sys_COMFF.X_DIMS_FREE,:, jj), ...
     UCOMFF(:,1:NUM_CTRL, jj), ...
     KCOMFF(:,sys_COMFF.X_DIMS_FREE,1:NUM_CTRL, jj), TraceCOMFF{jj,1}] = iLQG(dyn_COMFF, ...
                                                    XCOMFFinit(:,:, jj), ...
                                                    UCOMFFinit(:,1:NUM_CTRL, jj), Op);
    timeCOMFF(jj) = toc;
    XCOMFF(sys_COMFF.X_DIMS_FIXED,:, jj) = sys_COMFF.l_point(sys_COMFF.X_DIMS_FIXED) ...
                                              * ones(1, size(XCOMFF, 2));
    jj
end

% Torso second
A_ = [A(:,1:4) - B(:,1:2)*K_COMFF, A(:,5:6)];
B_ = [zeros(6,2), B(:,3:4)];
Q_ = sys.Q; Q_(1:4,1:4) = Q_(1:4,1:4) + K_COMFF'*sys.R(1:2,1:2)*K_COMFF;
R_ = sys.R;
[K_TorsoTS, S_TorsoTS, e_TorsoTS] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);
K_TorsoTS = K_TorsoTS(3:4,:);
K_COMFF = [K_COMFF, zeros(2,2)];

disp('T - Both')
sys_TorsoTS = sys;
sys_TorsoTS.U_DIMS_FREE = [3; 4];
sys_TorsoTS.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoTS.U_DIMS_FIXED(sys_TorsoTS.U_DIMS_FREE) = [];
sys_TorsoTS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_TorsoTS.X_DIMS_FIXED = [];

Op.lims = sys_TorsoTS.lims(sys_TorsoTS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [1;2;3;4];
XTorsoTS = zeros(length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoTS = zeros(length(sys_TorsoTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoTS = zeros(length(sys_TorsoTS.U_DIMS_FREE), length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XCOMFClose = nan(length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMFClose = zeros(length(sys_COMFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMFClose = zeros(length(sys_COMFF.U_DIMS_FREE), length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTorsoTSinit = nan(length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoTSinit = nan(length(sys_TorsoTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostTorsoTSinit = zeros(size(x_starts, 2), 1);
TraceTorsoTS = cell(size(x_starts, 2), 1);
timeTorsoTS = zeros(size(x_starts, 2), 1);

u0Free = sys_TorsoTS.u0(sys_TorsoTS.U_DIMS_FREE, 1);
u0Fixed = sys_TorsoTS.u0(sys_TorsoTS.U_DIMS_FIXED, 1);
goal = sys_TorsoTS.goal(sys_TorsoTS.X_DIMS_FREE, 1);
lims = sys_TorsoTS.lims(sys_TorsoTS.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_TorsoTS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_TorsoTS, x, u, k, K, xn, sys_TorsoTS.full_DDP);
    
    XTorsoTSinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UTorsoTSinit(:,ii, jj) = min(max(u0Free - K_TorsoTS*(XTorsoTSinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XTorsoTSinit(:,ii+1, jj) = f_Biped2DSecond_finite(sys_TorsoTS, XTorsoTSinit(:,ii, jj), UTorsoTSinit(:,ii, jj), u0Fixed, -K_COMFF, sys_COMFF.goal, sys_TorsoTS.dt);
        if (any(isnan(UTorsoTSinit(:,ii, jj))) || any(isnan(XTorsoTSinit(:,ii, jj))))
            disp('Check X, U TorsoTS');
        end
        CostTorsoTSinit(jj) = CostTorsoTSinit(jj) + discount*l_Biped2DSecond(sys_TorsoTS, XTorsoTSinit(:,ii, jj), UTorsoTSinit(:,ii, jj) - u0Free, zeros(length(u0Fixed),1), -K_COMFF, sys_COMFF.goal)*sys_TorsoTS.dt;
        discount = discount*sys_TorsoTS.gamma_;
    end
    CostTorsoTSinit(jj) = CostTorsoTSinit(jj) + discount*l_Biped2DSecond(sys_TorsoTS, XTorsoTSinit(:,NUM_CTRL+1, jj), zeros(length(u0Free),1), zeros(length(u0Fixed),1), -K_COMFF, sys_COMFF.goal)*sys_TorsoTS.dt;
    
    Op.cost = CostTorsoTSinit(jj);
    tic;
    [XTorsoTS(:,:, jj), ...
     UTorsoTS(:,1:NUM_CTRL, jj), ...
     KTorsoTS(:,:,1:NUM_CTRL, jj), ...
     UCOMFClose(:,:, jj), ...
     KCOMFClose(:,:,:, jj), ...
     XCOMFClose(:,:, jj), TraceTorsoTS{jj,1}] = iLQGSecondKDTree(dyn_TorsoTS, ...
                                        XTorsoTSinit(:,:, jj), ...
                                        UTorsoTSinit(:,1:NUM_CTRL, jj), ...
                                        UCOMFF(:,:, ceil(jj/size(theta_starts, 2))), ...
                                        KCOMFF(:,:,:, ceil(jj/size(theta_starts, 2))), ...
                                        XCOMFF(sys_COMFF.X_DIMS_FREE,:, ceil(jj/size(theta_starts, 2))), ...
                                        Op);
    timeTorsoTS(jj) = toc;
    jj

end

%% Torso - T, COM - F
disp('**** T - Torso, F - Both ****');
% Torso first
% A_ = A(5:6,5:6);
% B_ = B(5:6,3:4);

disp('T - Torso');
sys_TorsoTF = sys;
sys_TorsoTF.U_DIMS_FREE = [3; 4];
sys_TorsoTF.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoTF.U_DIMS_FIXED(sys_TorsoTF.U_DIMS_FREE) = [];
sys_TorsoTF.X_DIMS_FREE = [5; 6];
sys_TorsoTF.X_DIMS_FIXED = linspace(1,6,6)';
sys_TorsoTF.X_DIMS_FIXED(sys_TorsoTF.X_DIMS_FREE) = [];

u0 = sys.u0;
u0(sys_TorsoTF.U_DIMS_FIXED) = 0;
A_ = eval(subs(state_dyn_x, [x; u], [sys.l_point; u0]));
A_ = A_(5:6,5:6);
B_ = eval(subs(act_dyn_u, [x; u], [sys.l_point; u0]));
B_ = B_(5:6,3:4);
Q_ = sys.Q(5:6,5:6);
R_ = sys.R(3:4,3:4);
[K_TorsoTF, S_TorsoTF, e_TorsoTF] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);

Op.lims = sys_TorsoTF.lims(sys_TorsoTF.U_DIMS_FREE, :);
XTorsoTF = zeros(length(sys_TorsoTF.X_DIMS_FREE) + length(sys_TorsoTF.X_DIMS_FIXED), NUM_CTRL+1, size(theta_starts, 2));
UTorsoTF = zeros(length(sys_TorsoTF.U_DIMS_FREE), NUM_CTRL+1, size(theta_starts, 2));
KTorsoTF = zeros(length(sys_TorsoTF.U_DIMS_FREE), length(sys_TorsoTF.X_DIMS_FREE) + length(sys_TorsoTF.X_DIMS_FIXED), ...
                 NUM_CTRL+1, size(theta_starts, 2));
XTorsoTFinit = nan(length(sys_TorsoTF.X_DIMS_FREE), NUM_CTRL+1, size(theta_starts, 2));
UTorsoTFinit = nan(length(sys_TorsoTF.U_DIMS_FREE), NUM_CTRL+1, size(theta_starts, 2));
CostTorsoTFinit = zeros(size(theta_starts, 2), 1);
TraceTorsoTF = cell(size(theta_starts, 2), 1);
timeTorsoTF = zeros(size(theta_starts, 2), 1);

u0 = sys_TorsoTF.u0(sys_TorsoTF.U_DIMS_FREE, 1);
goal = sys_TorsoTF.goal(sys_TorsoTF.X_DIMS_FREE, 1);
lims = sys_TorsoTF.lims(sys_TorsoTF.U_DIMS_FREE, :);
for jj=1:1:size(theta_starts, 2)
    dyn_TorsoTF = @(x, u, i) biped2d_dyn_first_cst(sys_TorsoTF, x, u, sys_TorsoTF.full_DDP);
    
    XTorsoTFinit(:,1, jj) = theta_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UTorsoTFinit(:,ii, jj) = min(max(u0 - K_TorsoTF*(XTorsoTFinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XTorsoTFinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_TorsoTF, XTorsoTFinit(:,ii, jj), UTorsoTFinit(:,ii, jj), sys_TorsoTF.dt);
        if (any(isnan(UTorsoTFinit(:,ii, jj))) || any(isnan(XTorsoTFinit(:,ii, jj))))
            disp('Check X, U TorsoTF');
        end
        CostTorsoTFinit(jj) = CostTorsoTFinit(jj) + discount*l_Biped2DFirst(sys_TorsoTF, XTorsoTFinit(:,ii, jj), UTorsoTFinit(:,ii, jj) - u0)*sys_TorsoTF.dt;
        discount = discount*sys_TorsoTF.gamma_;
    end
    CostTorsoTFinit(jj) = CostTorsoTFinit(jj) + discount*l_Biped2DFirst(sys_TorsoTF, XTorsoTFinit(:,NUM_CTRL+1, jj), zeros(2,1))*sys_TorsoTF.dt;
    
    Op.cost = CostTorsoTFinit(jj);
    tic;
    [XTorsoTF(sys_TorsoTF.X_DIMS_FREE,:, jj), ...
     UTorsoTF(:,1:NUM_CTRL, jj), ...
     KTorsoTF(:,sys_TorsoTF.X_DIMS_FREE,1:NUM_CTRL, jj), TraceTorsoTF{jj,1}] = iLQG(dyn_TorsoTF, ...
                                                    XTorsoTFinit(:,:, jj), ...
                                                    UTorsoTFinit(:,1:NUM_CTRL, jj), Op);
    timeTorsoTF(jj) = toc;
    XTorsoTF(sys_TorsoTF.X_DIMS_FIXED,:, jj) = sys_TorsoTF.l_point(sys_TorsoTF.X_DIMS_FIXED) ...
                                               * ones(1, size(XTorsoTF, 2));
    jj
end

% COM second
A_ = [A(:,1:4), A(:,5:6) - B(:,3:4)*K_TorsoTF];
B_ = [B(:,1:2), zeros(6,2)];
Q_ = sys.Q; Q_(5:6,5:6) = Q_(5:6,5:6) + K_TorsoTF'*sys.R(3:4,3:4)*K_TorsoTF;
R_ = sys.R;
[K_COMFS, S_COMFS, e_COMFS] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);
K_COMFS = K_COMFS(1:2,:);
K_TorsoTF = [zeros(2,4), K_TorsoTF];

disp('F - Both');
sys_COMFS = sys;
sys_COMFS.U_DIMS_FREE = [1; 2];
sys_COMFS.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMFS.U_DIMS_FIXED(sys_COMFS.U_DIMS_FREE) = [];
sys_COMFS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMFS.X_DIMS_FIXED = [];

Op.lims = sys_COMFS.lims(sys_COMFS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [5;6];
XCOMFS = zeros(length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMFS = zeros(length(sys_COMFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMFS = zeros(length(sys_COMFS.U_DIMS_FREE), length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTorsoTClose = nan(length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoTClose = zeros(length(sys_TorsoTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoTClose = zeros(length(sys_TorsoTF.U_DIMS_FREE), length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XCOMFSinit = nan(length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMFSinit = nan(length(sys_COMFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostCOMFSinit = zeros(size(x_starts, 2), 1);
TraceCOMFS = cell(size(x_starts, 2), 1);
timeCOMFS = zeros(size(x_starts, 2), 1);

u0Free = sys_COMFS.u0(sys_COMFS.U_DIMS_FREE, 1);
u0Fixed = sys_COMFS.u0(sys_COMFS.U_DIMS_FIXED, 1);
goal = sys_COMFS.goal(sys_COMFS.X_DIMS_FREE, 1);
lims = sys_COMFS.lims(sys_COMFS.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_COMFS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_COMFS, x, u, k, K, xn, sys_COMFS.full_DDP);
    XCOMFSinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UCOMFSinit(:,ii, jj) = min(max(u0Free - K_COMFS*(XCOMFSinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XCOMFSinit(:,ii+1, jj) = f_Biped2DSecond_finite(sys_COMFS, XCOMFSinit(:,ii, jj), UCOMFSinit(:,ii, jj), u0Fixed, -K_TorsoTF, sys_TorsoTF.goal, sys_COMFS.dt);
        if (any(isnan(UCOMFSinit(:,ii, jj))) || any(isnan(XCOMFSinit(:,ii, jj))))
            disp('Check X, U COMFS');
        end
        CostCOMFSinit(jj) = CostCOMFSinit(jj) + discount*l_Biped2DSecond(sys_COMFS, XCOMFSinit(:,ii, jj), UCOMFSinit(:,ii, jj) - u0Free, zeros(length(u0Fixed),1), -K_TorsoTF, sys_TorsoTF.goal)*sys_COMFS.dt;
        discount = discount*sys_COMFS.gamma_;
    end
    CostCOMFSinit(jj) = CostCOMFSinit(jj) + discount*l_Biped2DSecond(sys_COMFS, XCOMFSinit(:,NUM_CTRL+1, jj), zeros(length(u0Free),1), zeros(length(u0Fixed),1), -K_TorsoTF, sys_TorsoTF.goal)*sys_COMFS.dt;
    
    Op.cost = CostCOMFSinit(jj);
    tic;
    [XCOMFS(:,:, jj), ...
     UCOMFS(:,1:NUM_CTRL, jj), ...
     KCOMFS(:,:,1:NUM_CTRL, jj), ...
     UTorsoTClose(:,:, jj), ...
     KTorsoTClose(:,:,:, jj), ...
     XTorsoTClose(:,:, jj), TraceCOMFS{jj,1}] = iLQGSecondKDTree(dyn_COMFS, ...
                                        XCOMFSinit(:,:, jj), ...
                                        UCOMFSinit(:,1:NUM_CTRL, jj), ...
                                        UTorsoTF(:,:, mod(jj-1, size(theta_starts,2))+1), ...
                                        KTorsoTF(:,:,:, mod(jj-1, size(theta_starts,2))+1), ...
                                        XTorsoTF(sys_TorsoTF.X_DIMS_FREE,:, mod(jj-1, size(theta_starts,2))+1), ...
                                        Op);
     timeCOMFS(jj) = toc;
     jj
end

%% COM - T, Torso - F
disp('**** T - COM, F - Both ****');
% COM first
% A_ = A(1:4,1:4);
% B_ = B(1:4,3:4);

disp('T - COM');
sys_COMTF = sys;
sys_COMTF.U_DIMS_FREE = [3; 4];
sys_COMTF.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMTF.U_DIMS_FIXED(sys_COMTF.U_DIMS_FREE) = [];
sys_COMTF.X_DIMS_FREE = [1;2;3;4];
sys_COMTF.X_DIMS_FIXED = linspace(1,6,6)';
sys_COMTF.X_DIMS_FIXED(sys_COMFF.X_DIMS_FREE) = [];

u0 = sys.u0;
u0(sys_COMTF.U_DIMS_FIXED) = 0;
A_ = eval(subs(state_dyn_x, [x; u], [sys.l_point; u0]));
A_ = A_(1:4,1:4);
B_ = eval(subs(act_dyn_u, [x; u], [sys.l_point; u0]));
B_ = B_(1:4,3:4);
Q_ = sys.Q(1:4,1:4);
R_ = sys.R(3:4,3:4);
[K_COMTF, S_COMTF, e_COMTF] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);

Op.lims = sys_COMTF.lims(sys_COMTF.U_DIMS_FREE, :);
XCOMTF = zeros(length(sys_COMTF.X_DIMS_FREE) + length(sys_COMTF.X_DIMS_FIXED), NUM_CTRL+1, size(com_starts, 2));
UCOMTF = zeros(length(sys_COMTF.U_DIMS_FREE), NUM_CTRL+1, size(com_starts, 2));
KCOMTF = zeros(length(sys_COMTF.U_DIMS_FREE), length(sys_COMTF.X_DIMS_FREE) + length(sys_COMTF.X_DIMS_FIXED), ...
                NUM_CTRL+1, size(com_starts, 2));
XCOMTFinit = nan(length(sys_COMFF.X_DIMS_FREE), NUM_CTRL+1, size(com_starts, 2));
UCOMTFinit = nan(length(sys_COMFF.U_DIMS_FREE), NUM_CTRL+1, size(com_starts, 2));
CostCOMTFinit = zeros(size(com_starts, 2), 1);
TraceCOMTF = cell(size(com_starts, 2), 1);
timeCOMTF = zeros(size(com_starts, 2), 1);

u0 = sys_COMTF.u0(sys_COMTF.U_DIMS_FREE, 1);
goal = sys_COMTF.goal(sys_COMTF.X_DIMS_FREE, 1);
lims = sys_COMTF.lims(sys_COMTF.U_DIMS_FREE, :);
for jj=1:1:size(com_starts, 2)
    dyn_COMTF = @(x, u, i) biped2d_dyn_first_cst(sys_COMTF, x, u, sys_COMTF.full_DDP);
    
    XCOMTFinit(:,1, jj) = com_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UCOMTFinit(:,ii, jj) = min(max(u0 - K_COMTF*(XCOMTFinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XCOMTFinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_COMTF, XCOMTFinit(:,ii, jj), UCOMTFinit(:,ii, jj), sys_COMTF.dt);
        if (any(isnan(UCOMTFinit(:,ii, jj))) || any(isnan(XCOMTFinit(:,ii, jj))))
            disp('Check X, U COMTF');
        end
        CostCOMTFinit(jj) = CostCOMTFinit(jj) + discount*l_Biped2DFirst(sys_COMTF, XCOMTFinit(:,ii, jj), UCOMTFinit(:,ii, jj) - u0)*sys_COMTF.dt;
        discount = discount*sys_COMTF.gamma_;
    end
    CostCOMTFinit(jj) = CostCOMTFinit(jj) + discount*l_Biped2DFirst(sys_COMTF, XCOMTFinit(:,NUM_CTRL+1, jj), zeros(2,1))*sys_COMTF.dt;
    
    Op.cost = CostCOMTFinit(jj);
    tic;
    [XCOMTF(sys_COMTF.X_DIMS_FREE,:, jj), ...
     UCOMTF(:,1:NUM_CTRL, jj), ...
     KCOMTF(:,sys_COMTF.X_DIMS_FREE,1:NUM_CTRL, jj), TraceCOMTF{jj,1}] = iLQG(dyn_COMTF, ...
                                                    XCOMTFinit(:,:, jj), ...
                                                    UCOMTFinit(:,1:NUM_CTRL, jj), Op);
    timeCOMTF(jj) = toc;
    XCOMTF(sys_COMTF.X_DIMS_FIXED,:, jj) = sys_COMTF.l_point(sys_COMTF.X_DIMS_FIXED)...
                                              * ones(1, size(XCOMTF, 2));
    jj
end

% Torso second
A_ = [A(:,1:4) - B(:,3:4)*K_COMTF, A(:,5:6)];
B_ = [B(:,1:2), zeros(6,2)];
Q_ = sys.Q; Q_(1:4,1:4) = Q_(1:4,1:4) + K_COMTF'*sys.R(3:4,3:4)*K_COMTF;
R_ = sys.R;
[K_TorsoFS, S_TorsoFS, e_TorsoFS] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);
K_TorsoFS = K_TorsoFS(1:2,:);
K_COMTF = [K_COMTF, zeros(2,2)];

disp('F - Both');
sys_TorsoFS = sys;
sys_TorsoFS.U_DIMS_FREE = [1; 2];
sys_TorsoFS.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoFS.U_DIMS_FIXED(sys_TorsoFS.U_DIMS_FREE) = [];
sys_TorsoFS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_TorsoFS.X_DIMS_FIXED = [];

Op.lims = sys_TorsoFS.lims(sys_TorsoFS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [1;2;3;4];
XTorsoFS = zeros(length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoFS = zeros(length(sys_TorsoFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoFS = zeros(length(sys_TorsoFS.U_DIMS_FREE), length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XCOMTClose = nan(length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMTClose = zeros(length(sys_COMTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMTClose = zeros(length(sys_COMTF.U_DIMS_FREE), length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTorsoFSinit = nan(length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoFSinit = nan(length(sys_TorsoFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostTorsoFSinit = zeros(size(x_starts, 2), 1);
TraceTorsoFS = cell(size(x_starts, 2), 1);
timeTorsoFS = zeros(size(x_starts, 2), 1);

u0Free = sys_TorsoFS.u0(sys_TorsoFS.U_DIMS_FREE, 1);
u0Fixed = sys_TorsoFS.u0(sys_TorsoFS.U_DIMS_FIXED, 1);
goal = sys_TorsoFS.goal(sys_TorsoFS.X_DIMS_FREE, 1);
lims = sys_TorsoFS.lims(sys_TorsoFS.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_TorsoFS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_TorsoFS, x, u, k, K, xn, sys_TorsoFS.full_DDP);
    
    XTorsoFSinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UTorsoFSinit(:,ii, jj) = min(max(u0Free - K_TorsoFS*(XTorsoFSinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XTorsoFSinit(:,ii+1, jj) = f_Biped2DSecond_finite(sys_TorsoFS, XTorsoFSinit(:,ii, jj), UTorsoFSinit(:,ii, jj), u0Fixed, -K_COMTF, sys_COMTF.goal, sys_TorsoFS.dt);
        if (any(isnan(UTorsoFSinit(:,ii, jj))) || any(isnan(XTorsoFSinit(:,ii, jj))))
            disp('Check X, U TorsoFS');
        end
        CostTorsoFSinit(jj) = CostTorsoFSinit(jj) + discount*l_Biped2DSecond(sys_TorsoFS, XTorsoFSinit(:,ii, jj), UTorsoFSinit(:,ii, jj) - u0Free, zeros(length(u0Fixed),1), -K_COMTF, sys_COMTF.goal)*sys_TorsoFS.dt;
        discount = discount*sys_TorsoFS.gamma_;
    end
    CostTorsoFSinit(jj) = CostTorsoFSinit(jj) + discount*l_Biped2DSecond(sys_TorsoFS, XTorsoFSinit(:,NUM_CTRL+1, jj), zeros(length(u0Free),1), zeros(length(u0Fixed),1), -K_COMTF, sys_COMTF.goal)*sys_TorsoFS.dt;
    
    Op.cost = CostTorsoFSinit(jj);
    tic;
    [XTorsoFS(:,:, jj), ...
     UTorsoFS(:,1:NUM_CTRL, jj), ...
     KTorsoFS(:,:,1:NUM_CTRL, jj), ...
     UCOMTClose(:,:, jj), ...
     KCOMTClose(:,:,:, jj), ...
     XCOMTClose(:,:, jj), TraceTorsoFS{jj,1}] = iLQGSecondKDTree(dyn_TorsoFS, ...
                                        XTorsoFSinit(:,:, jj), ...
                                        UTorsoFSinit(:,1:NUM_CTRL, jj), ...
                                        UCOMTF(:,:, ceil(jj/size(theta_starts, 2))), ...
                                        KCOMTF(:,:,:, ceil(jj/size(theta_starts, 2))), ...
                                        XCOMTF(sys_COMTF.X_DIMS_FREE,:, ceil(jj/size(theta_starts, 2))), ...
                                        Op);
    timeTorsoFS(jj) = toc;
    jj
end

%% Torso - F, COM - T
disp('**** F - Torso, T - Both ****');
% Torso first
% A_ = A(5:6,5:6);
% B_ = B(5:6,1:2);

disp('F - Torso');
sys_TorsoFF = sys;
sys_TorsoFF.U_DIMS_FREE = [1; 2];
sys_TorsoFF.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoFF.U_DIMS_FIXED(sys_TorsoFF.U_DIMS_FREE) = [];
sys_TorsoFF.X_DIMS_FREE = [5; 6];
sys_TorsoFF.X_DIMS_FIXED = linspace(1,6,6)';
sys_TorsoFF.X_DIMS_FIXED(sys_TorsoFF.X_DIMS_FREE) = [];

u0 = sys.u0;
u0(sys_TorsoFF.U_DIMS_FIXED) = 0;
A_ = eval(subs(state_dyn_x, [x; u], [sys.l_point; u0]));
A_ = A_(5:6,5:6);
B_ = eval(subs(act_dyn_u, [x; u], [sys.l_point; u0]));
B_ = B_(5:6,1:2);
Q_ = sys.Q(5:6,5:6);
R_ = sys.R(1:2,1:2);
[K_TorsoFF, S_TorsoFF, e_TorsoFF] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);

Op.lims = sys_TorsoFF.lims(sys_TorsoFF.U_DIMS_FREE, :);
XTorsoFF = zeros(length(sys_TorsoFF.X_DIMS_FREE) + length(sys_TorsoFF.X_DIMS_FIXED), NUM_CTRL+1, size(theta_starts, 2));
UTorsoFF = zeros(length(sys_TorsoFF.U_DIMS_FREE), NUM_CTRL+1, size(theta_starts, 2));
KTorsoFF = zeros(length(sys_TorsoFF.U_DIMS_FREE), length(sys_TorsoFF.X_DIMS_FREE) + length(sys_TorsoFF.X_DIMS_FIXED), ...
                NUM_CTRL+1, size(theta_starts, 2));
XTorsoFFinit = nan(length(sys_TorsoFF.X_DIMS_FREE), NUM_CTRL+1, size(theta_starts, 2));
UTorsoFFinit = nan(length(sys_TorsoFF.U_DIMS_FREE), NUM_CTRL+1, size(theta_starts, 2));
CostTorsoFFinit = zeros(size(theta_starts, 2), 1);
TraceTorsoFF = cell(size(theta_starts, 2), 1);
timeTorsoFF = zeros(size(theta_starts, 2), 1);

u0 = sys_TorsoFF.u0(sys_TorsoFF.U_DIMS_FREE, 1);
goal = sys_TorsoFF.goal(sys_TorsoFF.X_DIMS_FREE, 1);
lims = sys_TorsoFF.lims(sys_TorsoFF.U_DIMS_FREE, :);
for jj=1:1:size(theta_starts, 2)
    dyn_TorsoFF = @(x, u, i) biped2d_dyn_first_cst(sys_TorsoFF, x, u, sys_TorsoFF.full_DDP);
    
    XTorsoFFinit(:,1, jj) = theta_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UTorsoFFinit(:,ii, jj) = min(max(u0 - K_TorsoFF*(XTorsoFFinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));
        XTorsoFFinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_TorsoFF, XTorsoFFinit(:,ii, jj), UTorsoFFinit(:,ii, jj), sys_TorsoFF.dt);
        if (any(isnan(UTorsoFFinit(:,ii, jj))) || any(isnan(XTorsoFFinit(:,ii, jj))))
            disp('Check X, U TorsoFF');
        end
        CostTorsoFFinit(jj) = CostTorsoFFinit(jj) + discount*l_Biped2DFirst(sys_TorsoFF, XTorsoFFinit(:,ii, jj), UTorsoFFinit(:,ii, jj) - u0)*sys_TorsoFF.dt;
        discount = discount*sys_TorsoFF.gamma_;
    end
    CostTorsoFFinit(jj) = CostTorsoFFinit(jj) + discount*l_Biped2DFirst(sys_TorsoFF, XTorsoFFinit(:,NUM_CTRL+1, jj), zeros(2,1))*sys_TorsoFF.dt;
    
    Op.cost = CostTorsoFFinit(jj);tic;
    [XTorsoFF(sys_TorsoFF.X_DIMS_FREE,:, jj), ...
     UTorsoFF(:,1:NUM_CTRL, jj), ...
     KTorsoFF(:,sys_TorsoFF.X_DIMS_FREE,1:NUM_CTRL, jj), TraceTorsoFF{jj,1}] = iLQG(dyn_TorsoFF, ...
                                                    XTorsoFFinit(:,:, jj), ...
                                                    UTorsoFFinit(:,1:NUM_CTRL, jj), Op);
    timeTorsoFF(jj) = toc;
    XTorsoFF(sys_TorsoFF.X_DIMS_FIXED,:, jj) = sys_TorsoFF.l_point(sys_TorsoFF.X_DIMS_FIXED)...
                                              * ones(1, size(XTorsoFF, 2));
    jj
end

% COM second
A_ = [A(:,1:4), A(:,5:6) - B(:,1:2)*K_TorsoFF];
B_ = [zeros(6,2), B(:,3:4)];
Q_ = sys.Q; Q_(5:6,5:6) = Q_(5:6,5:6) + K_TorsoFF'*sys.R(1:2,1:2)*K_TorsoFF;
R_ = sys.R;
[K_COMTS, S_COMTS, e_COMTS] = lqr(A_ - lambda_/2*eye(size(A_,1)), B_, Q_, R_);
K_COMTS = K_COMTS(3:4,:);
K_TorsoFF = [zeros(2,4), K_TorsoFF];

disp('T - Both');
sys_COMTS = sys;
sys_COMTS.U_DIMS_FREE = [3; 4];
sys_COMTS.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMTS.U_DIMS_FIXED(sys_COMTS.U_DIMS_FREE) = [];
sys_COMTS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMTS.X_DIMS_FIXED = [];

Op.lims = sys_COMTS.lims(sys_COMTS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [5;6];
XCOMTS = zeros(length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMTS = zeros(length(sys_COMTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMTS = zeros(length(sys_COMTS.U_DIMS_FREE), length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTorsoFClose = nan(length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoFClose = zeros(length(sys_TorsoFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoFClose = zeros(length(sys_TorsoFF.U_DIMS_FREE), length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XCOMTSinit = nan(length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMTSinit = nan(length(sys_COMTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
CostCOMTSinit = zeros(size(x_starts, 2), 1);
TraceCOMTS = cell(size(x_starts, 2), 1);
timeCOMTS = zeros(size(x_starts, 2), 1);

u0Free = sys_COMTS.u0(sys_COMTS.U_DIMS_FREE, 1);
u0Fixed = sys_COMTS.u0(sys_COMTS.U_DIMS_FIXED, 1);
goal = sys_COMTS.goal(sys_COMTS.X_DIMS_FREE, 1);
lims = sys_COMTS.lims(sys_COMTS.U_DIMS_FREE, :);
for jj=1:1:size(x_starts, 2)
    dyn_COMTS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_COMTS, x, u, k, K, xn, sys_COMTS.full_DDP);
    
    XCOMTSinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        UCOMTSinit(:,ii, jj) = min(max(u0Free - K_COMTS*(XCOMTSinit(:,ii, jj) - goal), lims(:,1)), lims(:,2));        
        XCOMTSinit(:,ii+1, jj) = f_Biped2DSecond_finite(sys_COMTS, XCOMTSinit(:,ii, jj), UCOMTSinit(:,ii, jj), u0Fixed, -K_TorsoFF, sys_TorsoFF.goal, sys_COMTS.dt);
        if (any(isnan(UCOMTSinit(:,ii, jj))) || any(isnan(XCOMTSinit(:,ii, jj))))
            disp('Check X, U COMTS');
        end
        CostCOMTSinit(jj) = CostCOMTSinit(jj) + discount*l_Biped2DSecond(sys_COMTS, XCOMTSinit(:,ii, jj), UCOMTSinit(:,ii, jj) - u0Free, zeros(length(u0Fixed),1), -K_TorsoFF, sys_TorsoFF.goal)*sys_COMTS.dt;
        discount = discount*sys_COMTS.gamma_;
    end
    CostCOMTSinit(jj) = CostCOMTSinit(jj) + discount*l_Biped2DSecond(sys_COMTS, XCOMTSinit(:,NUM_CTRL+1, jj), zeros(length(u0Free),1), zeros(length(u0Fixed),1), -K_TorsoFF, sys_TorsoFF.goal)*sys_COMTS.dt;
    
    Op.cost = CostCOMTSinit(jj);
    tic;
    [XCOMTS(:,:, jj), ...
     UCOMTS(:,1:NUM_CTRL, jj), ...
     KCOMTS(:,:,1:NUM_CTRL, jj), ...
     UTorsoFClose(:,:, jj), ...
     KTorsoFClose(:,:,:, jj), ...
     XTorsoFClose(:,:, jj), TraceCOMTS{jj,1}] = iLQGSecondKDTree(dyn_COMTS, ...
                                        XCOMTSinit(:,:, jj), ...
                                        UCOMTSinit(:,1:NUM_CTRL, jj), ...
                                        UTorsoFF(:,:, mod(jj-1, size(theta_starts,2))+1), ...
                                        KTorsoFF(:,:,:, mod(jj-1, size(theta_starts,2))+1), ...
                                        XTorsoFF(sys_TorsoFF.X_DIMS_FREE,:, mod(jj-1, size(theta_starts,2))+1), ...
                                        Op);
    timeCOMTS(jj) = toc;
    jj
end

%% Decoupled
% COM - F, Torso - T
disp('**** F - COM, T - Torso ****');
sys_COMFTorsoTDec = sys;
sys_COMFTorsoTDec.U_DIMS_FREE = [1;2;3;4];
sys_COMFTorsoTDec.U_DIMS_FIXED = [];
sys_COMFTorsoTDec.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMFTorsoTDec.X_DIMS_FIXED = [];
XCOMFTorsoTDec = zeros(length(sys_COMFTorsoTDec.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMFTorsoTDec = zeros(length(sys_COMFTorsoTDec.U_DIMS_FREE), NUM_CTRL, size(x_starts, 2));
timeCOMFTorsoT = zeros(size(x_starts, 2), 1);
U_DIMS = cell(2,1);
U_DIMS{1,1} = [1;2];
U_DIMS{2,1} = [3;4];

for jj=1:1:size(x_starts, 2)
    dyn_COMFTorsoTDec = @(x, u) ...
               biped2d_dyn_first_cst(sys_COMFTorsoTDec, x, u, sys_COMFTorsoTDec.full_DDP);
    tic;
    [XCOMFTorsoTDec(:,:, jj), UCOMFTorsoTDec(:,:, jj)] = ForwardPassDec([UCOMFF(:,:, ceil(jj/size(theta_starts, 2))); UTorsoTF(:,:, mod(jj-1, size(theta_starts,2))+1)], ...
                                                                        [KCOMFF(:,:,:, ceil(jj/size(theta_starts, 2))); KTorsoTF(:,:,:, mod(jj-1, size(theta_starts,2))+1)], ...
                                                                        [XCOMFF(:,:, ceil(jj/size(theta_starts, 2))); XTorsoTF(:,:, mod(jj-1, size(theta_starts,2))+1)], ...
                                                                        U_DIMS, ...
                                                                         x_starts(:, jj), dyn_COMFTorsoTDec);
    timeCOMFTorsoT(jj) = toc;
    jj
end

% COM - T, Torso - F
disp('**** T - COM, F - Torso ****');
sys_COMTTorsoFDec = sys;
sys_COMTTorsoFDec.U_DIMS_FREE = [1;2;3;4];
sys_COMTTorsoFDec.U_DIMS_FIXED = [];
sys_COMTTorsoFDec.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMTTorsoFDec.X_DIMS_FIXED = [];
XCOMTTorsoFDec = zeros(length(sys_COMTTorsoFDec.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMTTorsoFDec = zeros(length(sys_COMTTorsoFDec.U_DIMS_FREE), NUM_CTRL, size(x_starts, 2));
timeCOMTTorsoF = zeros(size(x_starts, 2), 1);

for jj=1:1:size(x_starts, 2)
    dyn_COMTTorsoFDec = @(x, u) ...
               biped2d_dyn_first_cst(sys_COMTTorsoFDec, x, u, sys_COMTTorsoFDec.full_DDP);
    tic;
    [XCOMTTorsoFDec(:,:, jj), UCOMTTorsoFDec(:,:, jj)] = ForwardPassDec([UTorsoFF(:,:, mod(jj-1, size(theta_starts,2))+1); UCOMTF(:,:, ceil(jj/size(theta_starts, 2)))], ...
                                                                        [KTorsoFF(:,:,:, mod(jj-1, size(theta_starts,2))+1); KCOMTF(:,:,:, ceil(jj/size(theta_starts, 2)))], ...
                                                                        [XTorsoFF(:,:, mod(jj-1, size(theta_starts,2))+1); XCOMTF(:,:, ceil(jj/size(theta_starts, 2)))], ...
                                                                        U_DIMS, ...
                                                                        x_starts(:, jj), dyn_COMTTorsoFDec);
    timeCOMTTorsoF(jj) = toc;
    jj
end

save(strcat(save_dir, save_file, '.mat'));
