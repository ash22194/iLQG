function [x,u, u0] = test_cartpole
% A demo of iLQG/DDP for cartpole
clc; clear;
close all;
addpath('../iLQR_Matlab/cartPole/')
% Set full_DDP=true to compute 2nd order derivatives of the 
% dynamics. This will make iterations more expensive, but 
% final convergence will be much faster (quadratic)
full_DDP = false;
sys.mc = 5;
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81;
sys.dt = 0.001;
sys.gamma_ = 0.9995;
sys.Q = [25, 0, 0, 0;
         0, 0.02, 0, 0;
         0, 0, 25, 0;
         0, 0, 0, 0.02];
sys.R = [0.001, 0;
         0, 0.001];
sys.goal = [0; 0; pi; 0];

% set up the optimization problem
DYNCST  = @(x,u,i) cartpole_dyn_cst(sys, x, u, full_DDP);
global T;
T       = round(5 / sys.dt);              % horizon
global dt;
dt      = sys.dt;
global x0;  %[x,x-dot,theta,theta-dot]
x0      = [0; 0; 0; 0];   % initial state
global x_des;
x_des = sys.goal;

global u0; % initial controls
% TODO change this according to x0 and x_des?
u0      = zeros(2,T); % Just setting up shape here
u0(1,:) = 0.25*randn(1,T); % commanded force
u0(2,:) = 0.1*randn(1,T) + 3; % commanded torque
% u0(2,:) = 0.3*randn(1,T); % steering

Op.lims  = [-6 6;   
            -6 6];
Op.maxIter = 200;
Op.gamma_ = sys.gamma_;
global obs;
obs = [1; 0];   

% Initialize plot with start state, goal state, obstacles
% init_plot(x0,x_des,obs);

% Prepare trajectory visualization callback
% line_handle = line([0 0],[0 0],'color','b','linewidth',1.5);
% plotFn = @(x) traj_plot(x,line_handle);
% Op.plotFn = plotFn;

% === Run the optimization!
[x,u]= iLQG(DYNCST, x0, u0, Op);
plotForwardPass(x, [1,2;3,4], ["x", "x-dot"; "theta", "theta-dot"]);
% file_name = ['traj',datestr(now,'_mm-dd-yy_HH_MM')];
% save(['saved_trajectories/',file_name,'.mat'],'x','u','x0','x_des','dt','T');

end %test_cartpole