clear;
close all;
clc;

%% 
% Add path to relevant functions
addpath('../iLQR_Matlab/cartPole/');

sys.full_DDP = false;
sys.mc = 5;
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81;
sys.T = 5;
sys.dt = 0.001;
sys.gamma_ = 0.997;
sys.Q = [25, 0, 0, 0;
         0, 0.02, 0, 0;
         0, 0, 25, 0;
         0, 0, 0, 0.02];
sys.R = [0.005, 0;
         0, 0.005];
sys.goal = [0; 0; pi; 0];

% Dynamics and cost
sys.dynamics_discrete = @(X, U) f_CartPole_finite(sys, X, U, sys.dt);
% sys.dynamics_discrete = @(X, U) f_CartPole_WrapAround_finite(sys, X, U, sys.dt);
DYNCST  = @(x, u, i) cartpole_dyn_cst(sys, x, u, sys.full_DDP);

% Optimization parameters
Op.lims  = [-9 9;   
            -9 9];
Op.maxIter = 200;
Op.gamma_ = sys.gamma_;

% Linear system
x = sym('x', [4, 1]);
u = sym('u', [2, 1]);
zero_dyn = [x(2);
            (sys.mp*sys.l*(x(4)^2)*sin(x(3)) + 0.5*sys.mp*sys.g*sin(2*x(3)))/(sys.mc + sys.mp*sin(x(3))^2);
            x(4);
            (-0.5*sys.mp*sys.l*x(4)^2*sin(2*x(3)) - (sys.mc+sys.mp)*sys.g*sin(x(3)))/(sys.l*(sys.mc + sys.mp*sin(x(3))^2))];
act_dyn = [0;
           (u(1) - u(2)*cos(x(3))/sys.l)/(sys.mc + sys.mp*sin(x(3))^2);
           0;
           (u(2)/sys.l*(sys.mc/sys.mp+1) - u(1)*cos(x(3)))/(sys.l*(sys.mc + sys.mp*sin(x(3))^2))];
zero_dyn_x = jacobian(zero_dyn,x);
act_dyn_u = jacobian(act_dyn,u);
A = eval(subs(zero_dyn_x, x, sys.goal)); 
B = eval(subs(act_dyn_u, x, sys.goal));
lambda_ = (1 - sys.gamma_)/sys.dt;
[K_lqr, S_lqr, ~] = lqr(A-lambda_/2*eye(4), B, sys.Q, sys.R);

% Define starts
cart_starts = [-1, -0.5, 0, 0.5, 1;
                0,  0,   0, 0,   0];
pole_starts = [7*pi/4, 5*pi/4, 3*pi/4, pi/4, 0;
               0,      0,      0,      0,    0];
x_starts = nan(4, size(cart_starts,2)*size(pole_starts,2));
for ii=1:1:size(cart_starts, 2)
    for jj=1:1:size(pole_starts, 2)
        x_starts(:, (ii-1)*size(pole_starts, 2) + jj) = [cart_starts(:, ii); pole_starts(:, jj)];
    end
end
% x_starts = [-0.5;0;pi/4;0];

%% Run iLQG
NUM_CTRL = round(sys.T / sys.dt);
Uinit = nan(2, NUM_CTRL, size(x_starts, 2));
Xinit = nan(4, NUM_CTRL+1, size(x_starts, 2));
Costinit = nan(1, NUM_CTRL+1, size(x_starts, 2));
Ufinal = nan(2, NUM_CTRL, size(x_starts, 2));
Xfinal = nan(4, NUM_CTRL+1, size(x_starts, 2));
J_joint = nan(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    Xinit(:, 1, jj) = x_starts(:, jj);
    discount = 1;
    for ii = 1:1:NUM_CTRL
%         Uinit(:, ii, jj) = -K_lqr*(Xinit(:, ii, jj) - sys.goal);
        Uinit(:, ii, jj) = zeros(2,1);
        Costinit(:, ii, jj) = discount * l_CartPole(sys, Xinit(:, ii, jj), Uinit(:, ii, jj)) * sys.dt;
%         Costinit(:, ii, jj) = discount * l_CartPole_WrapAround(sys, Xinit(:, ii, jj), Uinit(:, ii, jj)) * sys.dt;
        discount = discount * sys.gamma_;
        Xinit(:, ii+1, jj) = sys.dynamics_discrete(Xinit(:, ii, jj), Uinit(:, ii, jj));
    end
    Costinit(:, NUM_CTRL+1, jj) = discount * l_CartPole(sys, Xinit(:, NUM_CTRL+1, jj), zeros(2,1)) * sys.dt;
%     Costinit(:, NUM_CTRL+1, jj) = discount * l_CartPole_WrapAround(sys, Xinit(:, NUM_CTRL+1, jj), zeros(2,1)) * sys.dt;
%     Op.cost = Costinit(:,:, jj);
    [Xfinal(:,:, jj), Ufinal(:,:, jj)] = iLQG(DYNCST, Xinit(:,1, jj), Uinit(:,:, jj), Op);
    J_joint(jj) = J_CartPole(sys, Xfinal(:,:, jj), Ufinal(:,:, jj));
%     J_joint(jj) = J_CartPole_WrapAround(sys, Xfinal(:,:, jj), Ufinal(:,:, jj));
    disp(strcat("Iteration : ", num2str(jj), ", J : ", num2str(J_joint(jj))));
end

% save(strcat('data/iLQG_mc=', num2str(sys.mc),',mp=',num2str(sys.mp),'_zeroinit.mat'));