clear;
close all;
clc;

%% 
% Add path to relevant functions
% addpath('../iLQR_Matlab/cartPole/');

sys.full_DDP = false;
sys.m = 5;
sys.T = 5;
sys.dt = 0.001;
sys.gamma_ = 0.9995;
sys.Q = [25, 0;
         0, 0.02];
sys.R = [0.005];
sys.goal = [0; 0];

% Dynamics and cost
DYNCST  = @(x, u, i) block_dyn_cst(sys, x, u, sys.full_DDP);

% Optimization parameters
Op.lims  = [-9 9];
Op.maxIter = 200;
Op.gamma_ = sys.gamma_;

% Define starts
x_starts = [-1, -0.5, 0, 0.5, 1;
             0,  0,   0, 0,   0];
x_starts = [-0.5;
             0];

%% Run iLQG
NUM_CTRL = round(sys.T / sys.dt);
Uinit = nan(1, NUM_CTRL, size(x_starts, 2));
Ufinal = nan(1, NUM_CTRL, size(x_starts, 2));
Xfinal = nan(2, NUM_CTRL+1, size(x_starts, 2));
J_joint = nan(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    Uinit(:,:,jj) = zeros(1, NUM_CTRL);
    [Xfinal(:,:, jj), Ufinal(:,:, jj)] = iLQG(DYNCST, x_starts(:, jj), Uinit(:,:, jj), Op);
    
end