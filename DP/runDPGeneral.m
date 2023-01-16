clear;
close all;
clc;

num_gpus = gpuDeviceCount();
gpu_id = 0;
max_avail_memory = 0;
for gg=1:1:num_gpus
    g = gpuDevice(gg);
    if (g.AvailableMemory > max_avail_memory)
        gpu_id = gg;
        max_avail_memory = g.AvailableMemory;
    end
end
g = gpuDevice(gpu_id);
reset(g);
fprintf('Using GPU : %d\n', gpu_id);

%% 

system_name = 'cartpole';
restoredefaultpath();
addpath('systems');
addpath(strcat('systems/', system_name));

decompositions = load(strcat('../iLQG_boxQP/data/multitraj/', system_name, '_pareto/summary.mat'));
sys = decompositions.sys;
sys.name = strcat(system_name, '_pareto');

Op.max_iter = sys.max_iter;
Op.max_policy_iter = sys.max_policy_iter;
Op.num_points = sys.num_points;
Op.num_action_samples = sys.num_action_samples;
Op.gtol = 1e-5;
Op.u_mean_tol = (sys.lims(:,2) - sys.lims(:,1)) * 2e-6;
Op.u_max_tol = (sys.lims(:,2) - sys.lims(:,1)) / 12;
Op.save_dir = 'data';
Op.reuse_policy = true;

if (size(decompositions.policy_decompositions, 1) > 44)
    policy_decompositions = decompositions.policy_decompositions_pareto_front;
else
    policy_decompositions = decompositions.policy_decompositions;
end

policies = cell(size(policy_decompositions, 1), 1);
value = cell(size(policy_decompositions, 1), 1);
info = cell(size(policy_decompositions, 1), 1);

for dd=1:1:size(policy_decompositions, 1)
    p = policy_decompositions{dd,1};
    s = policy_decompositions{dd,2};
    sys.decomposition_id = dd;
    fprintf('Decomposition : %d / %d\n', dd, size(policy_decompositions, 1));
    
    [policies{dd,1}, value{dd,1}, info{dd,1}] = dp_decomposition_gpu(sys, Op, p, s);
    info{dd,1} = rmfield(info{dd,1}, 'state_grid');
end

%% Joint 

p_joint = [zeros(sys.U_DIMS, 1), ones(sys.U_DIMS, 1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;

fprintf('Joint\n');
[policies_joint, value_joint, info_joint] = dp_decomposition_gpu(sys, Op, p_joint, s_joint);
info_joint = rmfield(info_joint, 'state_grid');

% Compute Error
% Find states within the bounds
grid_points = cell(sys.X_DIMS, 1);
for xx=1:1:sys.X_DIMS
    grid_points{xx} = linspace(sys.limits(xx,1), sys.limits(xx,2), sys.num_points(xx));
end
[grid_points{:}] = ndgrid(grid_points{:});

valid_points = true(size(grid_points{1}));
for xx=1:1:sys.X_DIMS
    valid_points = valid_points & ((grid_points{xx} <= sys.state_bounds(xx,2))...
                                   & (grid_points{xx} >= sys.state_bounds(xx,1)));
end

err_dp = cellfun(@(x) mean(x(valid_points) - value_joint(valid_points), 'all'), value);