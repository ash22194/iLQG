clear;
close all;
clc;

%%

restoredefaultpath;
system_name = 'quadcopter';
addpath('general');
addpath(strcat('new_systems/', system_name));
load(strcat('new_systems/', system_name, '/sys.mat'), 'sys');

sys.m = 0.5;
sys.I = diag([4.86*1e-3; 4.86*1e-3; 8.8*1e-3]);
sys.l = 0.225;
sys.g = 9.81;
sys.bk = 1.14*1e-7/(2.98*1e-6); % tau/f
sys.T = 6;
sys.dt = 0.001;
sys.X_DIMS = 10; % z, ro, pi, ya, vx, vy, vz, vr0, vpi, vya
sys.U_DIMS = 4;  % f1, f2, f3, f4 

sys.l_point = [1; zeros(9,1)];
sys.goal = sys.l_point;
sys.u0 = (sys.m*sys.g/sys.U_DIMS)*ones(sys.U_DIMS, 1);
sys.fxu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
sys.lims = [zeros(sys.U_DIMS, 1), 0.5*sys.m*sys.g*ones(sys.U_DIMS, 1)];

sys.gamma_ = 0.999;
sys.lambda_ = (1 - sys.gamma_)/sys.dt;
sys.Q = diag([1, 2, 2, 0.25, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01]);
sys.R = 0.002*eye(4);
sys.full_DDP = false;

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
% Op.Alpha = [1];

roll_pitch_starts = [-2*pi/8, 2*pi/8;
                      2*pi/8, 2*pi/8;
                      0,      0];
linear_velocity_starts = [0.5, 0.5, 0;
                          0.5, 0, 0.5];
angular_velocity_starts = [-0.5, 0.5, 0;
                            0.5, 0.5, 0;
                            0.5, 0, 0.5;
                            0, 0.5, 0.5];
height_starts = [0.6; 1.4];
starts = zeros(sys.X_DIMS, size(roll_pitch_starts, 1) ...
                           *size(linear_velocity_starts, 1) ...
                           *size(angular_velocity_starts, 1) ...
                           *size(height_starts, 1));
count = 0;
for rp=1:1:size(roll_pitch_starts,1)
    for lv=1:1:size(linear_velocity_starts, 1)
        for av=1:1:size(angular_velocity_starts, 1)
            for h=1:1:size(height_starts, 1)
                count = count + 1;
                starts(1, count) = height_starts(h, :)';
                starts(2:3, count) = roll_pitch_starts(rp, :)';
                starts(5:7, count) = linear_velocity_starts(lv, :)';
                starts(8:10, count) = angular_velocity_starts(av, :)';
            end
        end
    end
end

%% Test Decompositions

decompositions_file = strcat('../GA/data/', system_name, 'ParetoFront.mat');
load(decompositions_file, 'u_x');

Xd = zeros(sys.X_DIMS, round(sys.T / sys.dt)+1, size(starts, 2), size(u_x, 1));
Ud = zeros(sys.U_DIMS, round(sys.T / sys.dt), size(starts, 2), size(u_x, 1));
cd = zeros(1, size(starts, 2), size(u_x, 1));
for d=1:1:size(u_x, 1)
    
    fprintf('Decomposition : %d / %d\n', d, size(u_x, 1));
    
    p = reshape(u_x(d, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s = reshape(u_x(d, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    [Xd(:,:,:,d), Ud(:,:,:,d), cd(:,:,d)] = ilqg_decomposition(sys, Op, p, s, starts);
    
end

%% Test Joint Optimization

p_joint = [0, 1;
           0, 1;
           0, 1;
           0, 1];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);

[Xjoint, Ujoint, cjoint] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);
