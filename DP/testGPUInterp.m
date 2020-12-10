clear;
close all;
clc;

%% 

addpath('cuda');
num_dims = 8;
num_points = [9, 9, 9, 9, 7, 7, 7, 7]';
if (num_dims==1)
    values = rand(num_points, 1);
else
    values = rand(num_points');
end
grid_indices = cell(num_dims, 1);
grid = cell(num_dims, 1);
for xxi = 1:1:num_dims
    xx = num_points(xxi);
    grid_indices{xxi} = linspace(0, 1, num_points(xxi));
end
[grid{:}] = ndgrid(grid_indices{:});

num_iter = 5;
query = cell(num_iter, 1);
for jj=1:1:num_iter
    query{jj} = cell(num_dims, 1);
    for xxi=1:1:num_dims
        query{jj}{xxi} = rand(size(values));
    end
end

% Interpolation on CPU
tic;
G = griddedInterpolant(grid{:}, values);
for ii=1:1:num_iter
    grid_query_cpu = G(query{ii}{:});
end
time_cpu = toc;

grid = cellfun(@(x) gpuArray(x), grid, 'UniformOutput', false);
for jj=1:1:num_iter
    query{jj} = cellfun(@(x) gpuArray(x), query{jj}, 'UniformOutput', false);
end
values = gpuArray(values);

% Interpolation on GPU using in-built MATLAB function
tic;
for ii=1:1:num_iter
    if (num_dims > 5)
        grid_query_gpu = grid_query_cpu;
    else
        grid_query_gpu = interpn(grid{:}, values, query{ii}{:});
    end
end
time_gpu = toc;

% Interpolation when feval is called inside the function
tic;
for ii=1:1:num_iter
    grid_query_custom = interpn_gpukernel(grid{:}, values, query{ii}{:});
end
time_custom = toc;

% Interpolation when feval is called from the main script
nlinear = sprintf('calc_average%d',num_dims);
k = parallel.gpu.CUDAKernel(strcat('cuda/', nlinear, '.ptx'), ...
                            strcat('cuda/', nlinear, '.cu'), ...
                            nlinear);
k.ThreadBlockSize = 896;
k.GridSize = 1024;
Nx = size(grid{1}); % Number of grid points in each dimension
x1 = cellfun(@(y) y(1), grid); % Min value in each dimension
dx = cellfun(@(y) y(end) - y(1), grid) ./ (Nx(1:num_dims)' - 1); % Discretizatio in each dimension
corners = 0:1:(2^num_dims - 1);
cb = zeros(num_dims, 2^num_dims);
for ii=1:1:num_dims
    r = floor(corners / 2);
    cb(ii, :) = corners - 2*r;
    corners = r;
end
cb = num2cell(flip(cb + 1,1), 2);
corners_index = sub2ind(Nx, cb{:}) - 1;
x1 = num2cell(x1(:));
Nx = Nx(1:num_dims);
Nx = num2cell(Nx(:));
dx = num2cell(dx(:));

time_feval = zeros(num_iter, 1);
for ii=1:1:num_iter
    tic;
    grid_query_custom_kernel = feval(k, query{ii}{:}, values, Nx{:}, dx{:}, x1{:}, corners_index);
    time_feval(ii) = toc;
end
time_custom_kernel = sum(time_feval);

fprintf('Time (CPU) : %.4f\nTime (GPU built-in) : %.4f\nTime (GPU custom) : %.4f\nTime (GPU kernel) : %.4f\n', ...
         time_cpu, time_gpu, time_custom, time_custom_kernel);

% Calculate deviation in interpolation value on GPU vs CPU using different methods

max_discrepancy_gpu = max(abs(grid_query_gpu - grid_query_cpu), [], 'all');
max_discrepancy_custom = max(abs(grid_query_custom - grid_query_cpu), [], 'all');
max_discrepancy_custom_kernel = max(abs(grid_query_custom_kernel - grid_query_cpu), [], 'all');

fprintf('Discrepancy (GPU built-in) : %.5f\nDiscrepancy (GPU custom) : %.5f\nDiscrepancy (GPU kernel) : %.5f\n', ...
        max_discrepancy_gpu, max_discrepancy_custom, max_discrepancy_custom_kernel);