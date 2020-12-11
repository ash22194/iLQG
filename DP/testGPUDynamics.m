clear;
close all;
clc;
g = gpuDevice();
reset(g);

%%

system_name = 'manipulator4dof';
restoredefaultpath();
addpath(strcat('systems/', system_name));
addpath('systems');
load(strcat('data/', system_name, 'System.mat'));

% System parameters
sys.m = [5.4; 1.8; 0.6; 0.2]; % kg
sys.l = [0.2; 0.5; 0.25; 0.125]; % m
sys.limits = [pi/2, 3*pi/2; repmat([-pi/2, pi/2], 3, 1); repmat([-1.5, 1.5], 4, 1)];
sys.lims = [-24, 24; -15, 15; -7.5, 7.5; -1, 1]; % action limits
sys.g = 9.81; % m/s^2
M = num2cell(sys.m);
L = num2cell(sys.l);
sys.X_DIMS_FREE = linspace(1,8,8);

%% Generate grid of states
sys.num_points = [9,9,9,9,7,7,7,7];
grid_indices = cell(8, 1);
grid = cell(8, 1);
for xxi = 1:1:8
    grid_indices{xxi} = gpuArray(linspace(sys.limits(xxi,1), ...
                                          sys.limits(xxi,2), ...
                                          sys.num_points(xxi)));
end
[grid{:}] = ndgrid(grid_indices{:});

actions = cell(4, 1);
for uui=1:1:4
    actions{uui} = sys.lims(uui,1) + (sys.lims(uui,2) - sys.lims(uui,1)) ...
                                        * rand(size(grid{1}), 'gpuArray');
end

num_iter = 5;
time_arrayfun = 0;
time_kernel = 0;
time_kernel_wrapped = 0;
%% Dynamics on GPU using arrayfun
% disp('Running arrayfun...');
sys.kernel = parallel.gpu.CUDAKernel('systems/manipulator4dof/dyn.ptx', ...
                            'systems/manipulator4dof/dyn.cu', ...
                            'dyn');
sys.kernel.ThreadBlockSize = 256;
sys.kernel.GridSize = 1024;

% grid_arrayfun = dyn_finite_rk4_gpu(sys, grid, actions, 1e-3);
% f1 = @() dyn_finite_rk4_gpu(sys, grid, actions, 1e-3);
% times_arrayfun = zeros(num_iter, 1);
% for tt=1:1:num_iter
%     times_arrayfun(tt) = gputimeit(f1);
%     fprintf('Done iter : %d, Time : %.5f\n', tt, times_arrayfun(tt));
% end
% time_arrayfun = sum(times_arrayfun);

%% Dynamics on GPU using mex
% disp('Running mex...')
% profile on;
% grid_gpumex = dyn_finite_rk4_gpumex(sys, grid, actions, 1e-3);
% profile off;
% profsave(profile('info'), 'data/dyn_gpu_mex');
% f2 = @() dyn_finite_rk4_gpumex(sys, grid, actions, 1e-3);
% times_mex = zeros(num_iter, 1);
% for tt=1:1:num_iter
%     times_mex(tt) = gputimeit(f2);
%     fprintf('Done iter : %d, Time : %.5f\n', tt, times_mex(tt));
% end
% time_mex = sum(times_mex);

%% Dynamics on GPU using CUDA Kernel
disp('Running kernel...');
profile on;
grid_cudakernel = dyn_finite_rk4_gpukernel(sys, grid, actions, 1e-3);
profile off;
profsave(profile('info'), 'data/dyn_gpu_kernel');
f3 = @() dyn_finite_rk4_gpukernel(sys, grid, actions, 1e-3);
times_kernel = zeros(num_iter, 1);
for tt=1:1:num_iter
    times_kernel(tt) = gputimeit(f3);
    fprintf('Done iter : %d, Time : %.5f\n', tt, times_kernel(tt));
end
time_kernel = sum(times_kernel);

% fprintf('Time arrayfun : %.3f\nTime mex : %.5f\nTime kernel : %.5f\n', ...
%         time_arrayfun, time_mex, time_kernel);

%% Check discrepancy
% discrepancy_mex = zeros(8,1);
% discrepancy_kernel = zeros(8,1);
% for xxi=1:1:8
%     discrepancy_mex(xxi) = max(gather(grid_gpumex{xxi} - grid_arrayfun{xxi}), [], 'all');
%     discrepancy_kernel(xxi) = max(gather(grid_cudakernel{xxi} - grid_arrayfun{xxi}), [], 'all');
% end

% disp(strcat('Discrepancy mex :', sprintf('%.3f', discrepancy_mex)));
% disp(strcat('Discrepancy kernel :', sprintf('%.3f', discrepancy_kernel)));