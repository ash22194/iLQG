clear;
close all;
clc;
g = gpuDevice();
reset(g);

%%

system_name = 'manipulator3dof';
restoredefaultpath();
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/', system_name, 'System.mat'));
mexcuda(strcat('systems/', system_name, '/dyn_mex_finite.cu'), '-R2018a', '-output', strcat('systems/', system_name, '/dyn_mex_finite')); 

% System parameters
if (strcmp(system_name, 'manipulator4dof'))
    sys.m = [5.4; 1.8; 0.6; 0.2]; % kg
    sys.l = [0.2; 0.5; 0.25; 0.125]; % m
    sys.limits = [pi/2, 3*pi/2; repmat([-pi/2, pi/2], 3, 1); repmat([-1.5, 1.5], 4, 1)];
    sys.lims = [-24, 24; -15, 15; -7.5, 7.5; -1, 1]; % action limits
    sys.l_point = [pi; zeros(7,1)];
    sys.X_DIMS = 8;
    sys.X_DIMS_FREE = linspace(1,8,8)';
    sys.X_DIMS_FIXED = linspace(1,8,8)';
    sys.X_DIMS_FIXED(sys.X_DIMS_FREE) = [];

    sys.U_DIMS = 4;
    sys.U_DIMS_FREE = [1;2];
    sys.U_DIMS_CONTROLLED = [];
    sys.U_DIMS_FIXED = linspace(1,4,4)';
    sys.U_DIMS_FIXED(cat(1, sys.U_DIMS_FREE, sys.U_DIMS_CONTROLLED)) = [];
    
    sys.num_points = [9,9,9,9,7,7,7,7];
elseif (strcmp(system_name, 'manipulator3dof'))
    sys.m = [2.5; 0.5; 0.1] * 1.1; % kg
    sys.l = [0.5; 0.25; 0.125]; % m
    sys.limits = [0, 2*pi; repmat([-pi, pi], 2, 1); repmat([-3, 3], 3, 1)];
    sys.lims = [-16, 16; -7.5, 7.5; -1, 1]; % action limits
    sys.l_point = [pi; zeros(5,1)];
    sys.X_DIMS = 6;
    sys.X_DIMS_FREE = linspace(1,6,6)';
    sys.X_DIMS_FIXED = linspace(1,6,6)';
    sys.X_DIMS_FIXED(sys.X_DIMS_FREE) = [];

    sys.U_DIMS = 3;
    sys.U_DIMS_FREE = [1;2;3];
    sys.U_DIMS_CONTROLLED = [];
    sys.U_DIMS_FIXED = linspace(1,3,3)';
    sys.U_DIMS_FIXED(cat(1, sys.U_DIMS_FREE, sys.U_DIMS_CONTROLLED)) = [];
    
    sys.num_points = [17,17,17,13,13,13];
end

sys.g = 9.81; % m/s^2
M = num2cell(sys.m);
L = num2cell(sys.l);

%% Generate grid of states
sys.grid_size = ones(sys.X_DIMS,1);
sys.grid_size(sys.X_DIMS_FREE) = sys.num_points(sys.X_DIMS_FREE);
sys.grid_size = int32(sys.grid_size);

sys.active_actions = ones(sys.U_DIMS,1);
sys.active_actions(sys.U_DIMS_FIXED) = 0;
sys.active_actions = int32(sys.active_actions);

grid_indices = cell(length(sys.X_DIMS_FREE), 1);
grid = cell(sys.X_DIMS, 1);
for xxi = 1:1:length(sys.X_DIMS_FREE)
    xx = sys.X_DIMS_FREE(xxi);
    grid_indices{xxi} = linspace(sys.limits(xx,1), ...
                                 sys.limits(xx,2), ...
                                 sys.num_points(xx));
end
[grid{sys.X_DIMS_FREE}] = ndgrid(grid_indices{:});

for xxi=1:1:length(sys.X_DIMS_FIXED)
    xx = sys.X_DIMS_FIXED(xxi);
    grid{xx} = sys.l_point(xx);
end

actions = cell(sys.U_DIMS, 1);
for uui=1:1:sys.U_DIMS
    if (sys.active_actions(uui) == 1)
        actions{uui} = sys.lims(uui,1) + (sys.lims(uui,2) - sys.lims(uui,1)) ...
                                        * rand(size(grid{sys.X_DIMS_FREE(1)}));
    else
        actions{uui} = 0;
    end
end

num_iter = 5;
time_arrayfun = 0;
time_mex = 0;
time_kernel = 0;
%% Dynamics on CPU
disp('Running cpu...');
tic;
grid_cpu = dyn_finite_rk4(sys, grid, actions, 1e-3);
time_cpu = toc

grid = cellfun(@(x) gpuArray(x), grid, 'UniformOutput', false);
actions = cellfun(@(x) gpuArray(x), actions, 'UniformOutput', false);

%% Dynamics on GPU using arrayfun
% disp('Running arrayfun...');
% sys.kernel = parallel.gpu.CUDAKernel('systems/manipulator4dof/dyn.ptx', ...
%                             'systems/manipulator4dof/dyn.cu', ...
%                             'dyn');
% sys.kernel.ThreadBlockSize = 256;
% sys.kernel.GridSize = 1024;

% grid_arrayfun = dyn_finite_rk4_gpu(sys, grid, actions, 1e-3);
% f1 = @() dyn_finite_rk4_gpu(sys, grid, actions, 1e-3);
% times_arrayfun = zeros(num_iter, 1);
% for tt=1:1:num_iter
%     times_arrayfun(tt) = gputimeit(f1);
%     fprintf('Done iter : %d, Time : %.5f\n', tt, times_arrayfun(tt));
% end
% time_arrayfun = sum(times_arrayfun);

%% RK4 on GPU using mex
disp('Running RK4 mex...')
grid_rk4mex = dyn_finite_rk4_mex(sys, grid, actions, 1e-3);
f2 = @() dyn_finite_rk4_mex(sys, grid, actions, 1e-3);
times_rk4mex = zeros(num_iter, 1);
for tt=1:1:num_iter
    tic;
    grid_rk4mex = f2();
    times_rk4mex(tt) = toc;
    fprintf('Done iter : %d, Time : %.5f\n', tt, times_rk4mex(tt));
end
time_rk4mex = sum(times_rk4mex);

%% Dynamics on GPU using mex
% disp('Running mex...')
% grid_gpumex = dyn_finite_rk4_gpumex(sys, grid, actions, 1e-3);
% f3 = @() dyn_finite_rk4_gpumex(sys, grid, actions, 1e-3);
% times_mex = zeros(num_iter, 1);
% for tt=1:1:num_iter
%     times_mex(tt) = gputimeit(f3);
%     fprintf('Done iter : %d, Time : %.5f\n', tt, times_mex(tt));
% end
% time_mex = sum(times_mex);

%% Dynamics on GPU using CUDA Kernel
% disp('Running kernel...');
% profile on;
% grid_cudakernel = dyn_finite_rk4_gpukernel(sys, grid, actions, 1e-3);
% profile off;
% profsave(profile('info'), 'data/dyn_gpu_kernel');
% f4 = @() dyn_finite_rk4_gpukernel(sys, grid, actions, 1e-3);
% times_kernel = zeros(num_iter, 1);
% for tt=1:1:num_iter
%     times_kernel(tt) = gputimeit(f4);
%     fprintf('Done iter : %d, Time : %.5f\n', tt, times_kernel(tt));
% end
% time_kernel = sum(times_kernel);

% fprintf('Time arrayfun : %.3f\nTime mex : %.5f\nTime kernel : %.5f\n', ...
%         time_arrayfun, time_mex, time_kernel);

%% Check discrepancy
discrepancy_mex = zeros(sys.X_DIMS,1);
discrepancy_id = zeros(sys.X_DIMS,1);
% discrepancy_kernel = zeros(sys.X_DIMS,1);
for xxi=1:1:sys.X_DIMS
    discrepancy = abs(gather(grid_rk4mex{xxi}) - grid_cpu{xxi});
    [discrepancy_mex(xxi), discrepancy_id(xxi)] = max(discrepancy(:));
%     discrepancy_kernel(xxi) = max(gather(grid_cudakernel{xxi} - grid_arrayfun{xxi}), [], 'all');
end

disp(strcat('Discrepancy mex :', sprintf('%.5f ', discrepancy_mex)));
% disp(strcat('Discrepancy kernel :', sprintf('%.3f', discrepancy_kernel)));

init_state = cellfun(@(x) x(discrepancy_id(5)), grid);
init_action = cellfun(@(x) x(discrepancy_id(5)), actions);
next_state_cpu = cellfun(@(x) x(discrepancy_id(5)), grid_cpu);
next_state_gpu = cellfun(@(x) x(discrepancy_id(5)), grid_rk4mex);

disp(strcat('Discrepancy state : ', sprintf('%.3f ', init_state)));
disp(strcat('Discrepancy action : ', sprintf('%.3f ', init_action)));
disp(strcat('Next state (CPU) : ', sprintf('%.3f ', next_state_cpu)));
disp(strcat('Next state (GPU) : ', sprintf('%.3f ', next_state_gpu)));
