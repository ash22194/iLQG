function [value, info] = trajpolicy_evaluation_gpu(sys, Op, sub_policies)
%% Parameters
    X_DIMS            = sys.X_DIMS;
    X_DIMS_FREE       = sys.X_DIMS_FREE;
    X_DIMS_FIXED      = sys.X_DIMS_FIXED;
    U_DIMS            = sys.U_DIMS;
    U_DIMS_FREE       = sys.U_DIMS_FREE;
    U_DIMS_FIXED      = sys.U_DIMS_FIXED;
    U_DIMS_CONTROLLED = sys.U_DIMS_CONTROLLED;
    limits            = sys.limits;   % state limits
    time_limits       = sys.limits_t; % time length (sec) of trajectory to track
    lims              = sys.lims;     % action limits
    goal              = sys.goal;     % trajectory to track
    u0                = sys.u0;
    Q                 = sys.Q;
    R                 = sys.R;
    gamma_            = sys.gamma_;
    dt                = sys.dt;

    max_iter          = Op.max_iter;
    gtol              = Op.gtol;
    num_points        = Op.num_points;
    num_points_t      = size(goal,2);

%% Initialize 
    % Create state grids
    x = cell(X_DIMS+1, 1);
    grid_indices = cell(length(X_DIMS_FREE)+1, 1);
    for xxi = 1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        grid_indices{xxi} = gpuArray(linspace(limits(xx,1), limits(xx,2), num_points(xx)));
    end
    grid_indices{end} = gpuArray(linspace(time_limits(1), time_limits(2), num_points_t));
    [x{[X_DIMS_FREE; X_DIMS+1]}] = ndgrid(grid_indices{:});
    % Values for the fixed states? - Trajectory value at the time instant?
    for ff=1:1:length(X_DIMS_FIXED)
        ffi = X_DIMS_FIXED(ff);
        x{ffi} = gpuArray(repmat(reshape(goal(ffi,:), [ones(1,length(X_DIMS_FREE)), num_points_t]), ...
                                 [num_points(X_DIMS_FREE), 1]));
    end
    time_indices = repmat(reshape(linspace(1, num_points_t, num_points_t), ...
                                  [ones(1,length(X_DIMS_FREE)), num_points_t]), ...
                          [num_points(X_DIMS_FREE), 1]);
    
    % Initialize policy grids - fixed, free and controlled
    u = cell(U_DIMS, 1);
    u(U_DIMS_FIXED) = num2cell(zeros(length(U_DIMS_FIXED), 1));
    for uui = 1:1:length(U_DIMS_FREE)
        uu = U_DIMS_FREE(uui);
        u{uu} = lims(uu, 1) + (lims(uu, 2) - lims(uu, 1)) ...
                              * rand([num_points(X_DIMS_FREE), num_points_t], 'gpuArray');
    end
    for jj=1:1:size(sub_policies, 1)
        U_SUBDIM = sub_policies{jj,1};
        X_SUBDIM = sub_policies{jj,2};
        [X_SUBDIM, ~] = find(X_DIMS_FREE == X_SUBDIM');
        X_SUBDIM_BAR = 1:length(X_DIMS_FREE);
        X_SUBDIM_BAR(X_SUBDIM) = [];
        
        subpolicy_size = [num_points(X_DIMS_FREE), num_points_t];
        subpolicy_size(X_SUBDIM_BAR) = 1;
        subpolicy_newsize = num_points(X_DIMS_FREE);
        subpolicy_newsize(X_SUBDIM) = 1;
        subpolicy_newsize = cat(2, subpolicy_newsize, 1);
        
        u(U_SUBDIM) = cellfun(@(x) gpuArray(repmat(reshape(x, subpolicy_size), subpolicy_newsize)), sub_policies{jj,3}, 'UniformOutput', false);
    end
    
    % Initialize cost functions
    % State cost assuming diagonal Q matrix
    cost_state = 0;
    for xxi = 1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        goal_xx = goal(xx, :);
        cost_state = cost_state + dt * Q(xx,xx) * (x{xx} - goal_xx(time_indices)).^2;
    end
    
    % Controlled action cost assuming diagonal R matrix
    cost_action_controlled = 0;
    for uui=1:1:length(U_DIMS_CONTROLLED)
        uu = U_DIMS_CONTROLLED(uui);
        cost_action_controlled = cost_action_controlled + dt * R(uu, uu) * (u{uu} - u0(uu)).^2;
    end
    cost_total = cost_state + cost_action_controlled;
    
    goal_grid = cell(length(X_DIMS_FREE), 1);
    for xxi=1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        [~, goal_grid{xxi}] = min(abs(grid_indices{xxi}' - goal(xx,:)), [], 1);
    end
    goal_grid = cat(1, goal_grid, {linspace(1, num_points_t, num_points_t)});
    goal_grid = sub2ind([num_points(X_DIMS_FREE),num_points_t], goal_grid{:});
    
    G_ = zeros(num_points(X_DIMS_FREE));
    G = ones(size(G_));
    [kernel_interp, kernel_inputs] = generate_interp_kernel(x{[X_DIMS_FREE; X_DIMS+1]});
    policy_iter = 0;
    time_total = 0;
    
    time_start = tic;
    x_ = dyn_finite_rk4(sys, x, u, dt); % Next state after a step
    while ((max(abs(G_ - G), [], 'all') > gtol) && (policy_iter < max_iter))
        % Iterate to estimate value function
        tic;
        policy_iter = policy_iter + 1;

        G = G_;
        Gnext = feval(kernel_interp, x_{[X_DIMS_FREE; X_DIMS+1]}, G_, kernel_inputs{:});
        Gnext(goal_grid) = 0;

        % Update value function
        G_ = cost_total + gamma_*Gnext;
        G_(goal_grid) = 0;
        
        time_iter = toc;
%         disp(strcat('Policy evaluation iter :', num2str(policy_iter), ', time : ', num2str(time_iter)));
    end
    time_total = toc(time_start);
    
    value = G_;
    info.state_grid = x([X_DIMS_FREE; X_DIMS+1]);
    info.policy_iter = policy_iter;
    info.time_total = time_total;
end

function [k, inputs] = generate_interp_kernel(varargin)
% x1, ..., xn generated from ndgrid

n = nargin;
grid_sizes = cell2mat(cellfun(@(x) size(x), varargin(:), 'UniformOutput', false));
size_compare = ones(n ,1) * grid_sizes(1,:) == grid_sizes;
assert(all(size_compare, 'all'), "All xi must have the same dimensions");

Nx = grid_sizes(1,:)';
x1 = cellfun(@(y) y(1), varargin(:)); % Min value in each dimension
dx = cellfun(@(y) y(end) - y(1), varargin(:)) ./ (Nx(1:n) - 1); % Discretizatio in each dimension

% Generate corners of the n dimensional hyper-cube
corners = 0:1:(2^n - 1);
cb = zeros(n, 2^n);
for ii=1:1:n
    r = floor(corners / 2);
    cb(ii, :) = corners - 2*r;
    corners = r;
end
cb = num2cell(flip(cb + 1, 1), 2);
corners_index = sub2ind(Nx, cb{:}) - 1;

grid_size = prod(Nx(1:n));
Nx = num2cell(Nx(1:n));
x1 = num2cell(x1);
dx = num2cell(dx);
inputs = cat(1, {grid_size}, Nx, dx, x1, {corners_index});

nlinear = sprintf('calc_average%d',n);
k = parallel.gpu.CUDAKernel(strcat('cuda/', nlinear, '.ptx'), ...
                            strcat('cuda/', nlinear, '.cu'), ...
                            nlinear);
k.ThreadBlockSize = 896;
k.GridSize = 1024;
end