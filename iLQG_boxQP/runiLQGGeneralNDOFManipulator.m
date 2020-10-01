clear;
close all;
clc;

%%

restoredefaultpath;
n = 2;
system_name = sprintf('manipulator%ddof', n);
addpath('general');
addpath(strcat('new_systems/', system_name));
load(strcat('new_systems/', system_name, '/sys.mat'), 'sys');

if (n==2)
    sys.m = [2.5; 0.5]; % kg
    sys.l = [0.5; 0.25]; % m
elseif (n==3)
    sys.m = [12.5; 1.25; 0.125]; % kg
    sys.l = [0.6;  0.2;  0.067]; % m
end
sys.g = 9.81; % m/s^2
sys.T = 4;
sys.dt = 0.001;
sys.X_DIMS = 2*sys.n; % [thi, ... dthi, ...]
sys.U_DIMS = sys.n;   % [taui] 

sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(1) = pi;
sys.goal = sys.l_point;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.fxfu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
sys.lims = inf*[-ones(sys.n,1), ones(sys.n,1)]; % action limits

sys.gamma_ = 0.997;
sys.lambda_ = (1 - sys.gamma_)/sys.dt;
sys.Q = diag([25*ones(1, sys.n), 0.02*ones(1, sys.n)]);
sys.R = diag(0.001*ones(1, sys.n));
sys.full_DDP = false;
sys.u0init = false;

% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
% Op.Alpha = [1];

theta_starts(:,:,1) = [2*pi/3, 2*pi/3, 4*pi/3, 4*pi/3;
                       -1,     1,      -1,     1];
theta_starts(:,:,2:n) = repmat([-pi/3, -pi/3, pi/3, pi/3;
                                -1,    1,     -1,   1], [1,1,n-1]);

starts = zeros(sys.X_DIMS, size(theta_starts, 2)^n);
for count = 1:1:size(starts,2)
    count_ = count;
    for jj=1:1:n
        index = 1 + mod(count_-1, size(theta_starts, 2));
        starts([jj; n+jj], count) = theta_starts(:, index, jj);
        count_ = (count_ - index)/size(theta_starts, 2);
    end
end

%% Test Decompositions

% p = [linspace(0,n-1,n)', ones(n,1)];
% s = repmat(eye(n), [1, n]);
% u_x = [reshape(p, 1, 2*sys.U_DIMS), reshape(s, 1, sys.U_DIMS*sys.X_DIMS)];
% 
% Xd = zeros(sys.X_DIMS, round(sys.T / sys.dt)+1, size(starts, 2), size(u_x, 1));
% Ud = zeros(sys.U_DIMS, round(sys.T / sys.dt), size(starts, 2), size(u_x, 1));
% cd = zeros(1, size(starts, 2), size(u_x, 1));
% for d=1:1:size(u_x, 1)
%     
%     fprintf('Decomposition : %d / %d\n', d, size(u_x, 1));
%     sys.decomposition_id = d;
%     p = reshape(u_x(d, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
%     s = reshape(u_x(d, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
%     [Xd(:,:,:,d), Ud(:,:,:,d), cd(:,:,d)] = ilqg_decomposition(sys, Op, p, s, starts);
%     
% end

%% Test Joint Optimization

p_joint = [zeros(n,1), ones(n,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjoint, Ujoint, cjoint] = ilqg_decomposition(sys, Op, p_joint, s_joint, starts);
