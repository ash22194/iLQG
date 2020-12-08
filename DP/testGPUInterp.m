clear;
close all;
clc;

%% 

addpath('cuda');
num_dims = 1;
num_points = [69]';
if (num_dims==1)
    blah = rand(num_points, 1);
else
    blah = rand(num_points');
end
grid_indices = cell(num_dims, 1);
grid = cell(num_dims, 1);
query = cell(num_dims, 1);
for xxi = 1:1:num_dims
    xx = num_points(xxi);
    grid_indices{xxi} = linspace(0, 1, num_points(xxi));
    query{xxi} = rand(size(blah));
end
[grid{:}] = ndgrid(grid_indices{:});

tic;
G = griddedInterpolant(grid{:}, blah);
for ii=1:1:100
    grid_query_interpn = G(query{:});
end
time_interpn = toc;

grid = cellfun(@(x) gpuArray(x), grid, 'UniformOutput', false);
query = cellfun(@(x) gpuArray(x), query, 'UniformOutput', false);
blah = gpuArray(blah);

tic;
for ii=1:1:100
    grid_query_custom_old = interpn(grid{:}, blah, query{:});
end
time_custom_old = toc;

nlinear = sprintf('calc_average%d',num_dims);
k = parallel.gpu.CUDAKernel(strcat('cuda/', nlinear, '.ptx'), ...
                            strcat('cuda/', nlinear, '.cu'), ...
                            nlinear);
k.ThreadBlockSize = 1024;
k.GridSize = 1024;
tic;
for ii=1:1:100
    grid_query_custom = interpn_gpukernel(grid{:}, blah, query{:}, k);
end
time_custom = toc;

max_discrepancy = max(abs(grid_query_custom - grid_query_interpn), [], 'all');