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
[K_joint, S_joint, e_joint] = lqr(A - lambda_/2*eye(size(A,1)), B, ...
                                  sys.Q, sys.R, ...
                                  zeros(size(A,1), size(B,2)));

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
% Op.Alpha = [1];

% Define starts
% com_pos = [0.85, 0.85, 0.9, 0.9, 0.95, 0.95;
%             0.1,  0.4, 0.1, 0.4,  0.1,  0.4];
com_pos = [0.92, 0.92, 1.0, 1.0;
           0.4,  0.3, 0.4, 0.3];
com_pos(2,:) = pi/2 + com_pos(2,:);
com_vel = [ 0.1, -0.1, 0.1, -0.1;
           -0.3, -0.3, -0.4, -0.4];
theta_starts = [-0.4,  -0.25, 0.25, 0.4;
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
Xinit = nan(6, NUM_CTRL+1, size(x_starts, 2));
Uinit = nan(4, NUM_CTRL+1, size(x_starts, 2));
Costinit = zeros(size(x_starts, 2), 1);
TraceJoint = cell(size(x_starts, 2), 1);
timeJoint = zeros(size(x_starts, 2), 1);

for jj=1:1:size(x_starts, 2)
    dyn_joint = @(x, u, i) biped2d_dyn_first_cst(sys_joint, x, u, sys_joint.full_DDP);
    
%     l2init = sqrt((sys.df + x_starts(1)*cos(x_starts(2)))^2 + (x_starts(1)*sin(x_starts(2)))^2);
%     alpha2init = acos((sys.df + x_starts(1)*cos(x_starts(2)))/l2init);
%     uinit = [sys.m*sys.g*cos(alpha2init)/sin(x_starts(2) - alpha2init);
%              -sys.m*sys.g*cos(x_starts(2))/sin(x_starts(2) - alpha2init);
%               0;
%               0];
%     Uinit = repmat(uinit, [1, NUM_CTRL]);
    Xinit(:,1, jj) = x_starts(:,jj);
    discount = 1;
    for ii=1:1:NUM_CTRL
        Uinit(:,ii, jj) = min(max(sys_joint.u0 - K_joint*(Xinit(:,ii, jj) - sys_joint.goal), sys_joint.lims(:,1)), sys_joint.lims(:,2));
        Xinit(:,ii+1, jj) = f_Biped2DFirst_finite(sys_joint, Xinit(:,ii, jj), Uinit(:,ii, jj), sys_joint.dt);
        Costinit(jj) = Costinit(jj) + discount*l_Biped2DFirst(sys_joint, Xinit(:,ii, jj), Uinit(:,ii, jj))*sys_joint.dt;
        discount = discount*sys_joint.gamma_;
    end
    Costinit(jj) = Costinit(jj) + discount*l_Biped2DFirst(sys_joint, Xinit(:,NUM_CTRL+1, jj), zeros(4,1))*sys_joint.dt;
%     figure;
%     subplot(4,1,1);
%     x_coord = Xinit(1,:).*cos(Xinit(2,:));
%     z_coord = Xinit(1,:).*sin(Xinit(2,:));
%     plot(x_coord, z_coord); xlabel('x'); ylabel('z');
%     subplot(4,1,2);
%     plot(x_coord, Xinit(3,:)); xlabel('x'); ylabel('x-dot');
%     subplot(4,1,3);
%     plot(z_coord, Xinit(4,:)); xlabel('z'); ylabel('z-dot');
%     subplot(4,1,4);
%     plot(Xinit(5,:), Xinit(6,:)); xlabel('theta'); ylabel('theta-dot');
    
%     [~, Xinit] = ode45(@(t,x) f_Biped2DFirst(sys_joint, x, sys.u0 - K_joint*(x - sys.goal)), ...
%                        linspace(0, sys.T, NUM_CTRL+1), x_starts(:,jj));
%     Xinit = Xinit';
%     Uinit = sys.u0 - K_joint*(Xinit - sys.goal);
%     discount = sys.gamma_.^linspace(0, NUM_CTRL, NUM_CTRL+1);
%     Costinit = discount.*(diag(Xinit'*sys.Q*Xinit) + diag(Uinit'*sys.R*Uinit));
    
    Op.cost = Costinit;
    tic;
    [XJoint(sys_joint.X_DIMS_FREE,:, jj), ...
     UJoint(:,1:NUM_CTRL, jj), ...
     KJoint(:,sys_joint.X_DIMS_FREE,1:NUM_CTRL, jj), TraceJoint{jj,1}] = iLQG(dyn_joint, ...
                                                    Xinit(:,:,jj), ... % x_starts(sys_joint.X_DIMS_FREE, jj), ...
                                                    Uinit(:,1:NUM_CTRL,jj), Op);
    timeJoint(jj) = toc;
    XJoint(sys_joint.X_DIMS_FIXED,:, jj) = repmat(sys_joint.l_point(sys_joint.X_DIMS_FIXED), [1, size(XJoint,2)]);
    jj
    XJoint(:,end,jj)
end
