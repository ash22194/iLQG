function [policy, info] = get_trajdp_solution_gpu(sys, Op, sub_policies)
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
    
    max_iter           = Op.max_iter;
    max_policy_iter    = Op.max_policy_iter;
    gtol               = Op.gtol;
    u_max_tol          = Op.u_max_tol;
    u_mean_tol         = Op.u_mean_tol;
    num_points         = Op.num_points;
    num_points_t       = size(goal,2);
    num_action_samples = prod(Op.num_action_samples(U_DIMS_FREE));
    
    % Logging
    save_dir       = Op.save_dir;
    subsystem_id   = strcat('U', sprintf('%d', sys.U_DIMS_FREE), '_X', sprintf('%d', sys.X_DIMS_FREE));
    save_file      = strcat(save_dir, '/', subsystem_id, '_interm.mat');
    if (isfield(Op, 'save_every')) 
        save_every = Op.save_every; 
    else 
        save_every = max_iter / 10; 
    end
    reuse_policy   = isfield(Op, 'reuse_policy') && (Op.reuse_policy) && (isfile(save_file));
    
%% Initialize 
    % Create state grids
    % Initialize policy grids - fixed, free and controlled
    x = cell(X_DIMS+1, 1);
    u = cell(U_DIMS, 1);
    u(U_DIMS_FIXED) = num2cell(gpuArray(zeros(length(U_DIMS_FIXED), 1)));
    grid_indices = cell(length(X_DIMS_FREE)+1, 1);
    for xxi = 1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        grid_indices{xxi} = gpuArray(linspace(limits(xx,1), limits(xx,2), num_points(xx)));
    end
    grid_indices{end} = gpuArray(linspace(time_limits(1), time_limits(2), num_points_t));

    if (reuse_policy)
        decomposition  = load(save_file);
        u(U_DIMS_FREE) = decomposition.policy;
        info           = decomposition.info;
        x(X_DIMS_FREE) = info.state_grid(1:length(X_DIMS_FREE));
        x(end)         = info.state_grid(end);
        G_             = info.value;
    else
        for uui = 1:1:length(U_DIMS_FREE)
            uu = U_DIMS_FREE(uui);
            u{uu} = lims(uu, 1) + (lims(uu, 2) - lims(uu, 1)) ...
                                  * rand([num_points(X_DIMS_FREE), num_points_t], ...
                                         'gpuArray');
        end
        [x{[X_DIMS_FREE; X_DIMS+1]}] = ndgrid(grid_indices{:});
        G_ = gpuArray(zeros([num_points(X_DIMS_FREE), num_points_t]));

        info.iter = 0;
        info.time_total = 0;
        info.time_policy_eval = 0;
        info.time_policy_update = 0;
    end
    % Values for the fixed states? - Trajectory value at the time instant?
    for ff=1:1:length(X_DIMS_FIXED)
        ffi = X_DIMS_FIXED(ff);
        x{ffi} = gpuArray(repmat(reshape(goal(ffi,:), [ones(1,length(X_DIMS_FREE)), num_points_t]), ...
                                [num_points(X_DIMS_FREE), 1]));
    end
    time_indices = repmat(reshape(linspace(1, num_points_t, num_points_t), ...
                                  [ones(1,length(X_DIMS_FREE)), num_points_t]), ...
                          [num_points(X_DIMS_FREE), 1]);
    
    % Resize sub_policies
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
        
        u(U_SUBDIM) = cellfun(@(x) gpuArray(repmat(reshape(x, subpolicy_size), subpolicy_newsize)), ...
                              sub_policies{jj,3}, 'UniformOutput', false);
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
    cost_fixed = cost_state + cost_action_controlled;
    
    goal_grid = cell(length(X_DIMS_FREE), 1);
    for xxi=1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        [~, goal_grid{xxi}] = min(abs(grid_indices{xxi}' - goal(xx,:)), [], 1);
    end
    goal_grid = cat(1, goal_grid, {linspace(1, num_points_t, num_points_t)});
    goal_grid = sub2ind(size(G_), goal_grid{:});
    assert(numel(goal_grid)==num_points_t, 'Check number of goal states!');
    
    [kernel_interp, kernel_inputs] = generate_interp_kernel(x{[X_DIMS_FREE; X_DIMS+1]});
    constant_policy_count = 0;
    iter = info.iter;
    time_total = info.time_total;
    time_policy_eval = info.time_policy_eval;
    time_policy_update = info.time_policy_update;
    
%% Policy Iteration
    time_start = tic;
    while (iter < max_iter)
        iter = iter + 1;
        if (mod(iter, save_every)==0)
            policy = u(U_DIMS_FREE);
            info.value = G_;
            info.state_grid = x([X_DIMS_FREE; X_DIMS+1]);
            info.iter = iter;
            info.time_total = time_total;
            info.time_policy_eval = time_policy_eval;
            info.time_policy_update = time_policy_update;
            sys_ = sys;
            disp(strcat('Saving at iter', num2str(iter)));
            save(save_file, 'policy', 'info', 'sys_', '-v7.3', '-nocompression');
        end
        
        G = gpuArray(ones(size(G_)));
        policy_iter = 0;
        cost_action = 0;
        x_ = dyn_finite_rk4(sys, x, u, dt);
        for uui = 1:1:length(U_DIMS_FREE)
            uu = U_DIMS_FREE(uui);
            cost_action = cost_action + dt * R(uu, uu) * (u{uu} - u0(uu)).^2;
        end
        
        cost_total = cost_fixed + cost_action;
        time_policy_eval_start = tic;
        % Compute value function for current policy
         while ((max(abs(G_ - G), [], 'all') > gtol) && (policy_iter < max_policy_iter))
            % Iterate to estimate value function
            G = G_;
            policy_iter = policy_iter + 1;
            Gnext = feval(kernel_interp, x_{[X_DIMS_FREE; X_DIMS+1]}, G_, kernel_inputs{:});
            Gnext(goal_grid) = 0;
            
            % Update value function
            G_ = cost_total + gamma_*Gnext;
            G_(goal_grid) = 0;
        end
        time_policy_eval = time_policy_eval + toc(time_policy_eval_start);
        clear 'G';

        % Update policy
        % Compute Q(x,u) values g(x,u) + G(x+udt) for different actions when value
        % function G has converged and greedily update the policy i.e. choose
        % action u* at state x such that u* = argmin Q(x,u)
        minG = G_;
        unew = u;
        u_ = u;
        time_policy_update_start = tic;
        for i = 1:1:num_action_samples
            cost_action_ = 0;
            for uui=1:1:length(U_DIMS_FREE)
                uu = U_DIMS_FREE(uui);
                u_{uu} = lims(uu, 1) + (lims(uu, 2) - lims(uu, 1)) * rand(size(G_), 'gpuArray');
                cost_action_ = cost_action_ + dt * R(uu, uu) * (u_{uu} - u0(uu)).^2;
            end
            
            x__ = dyn_finite_rk4(sys, x, u_, dt);
            
            Gnext_ = feval(kernel_interp, x__{[X_DIMS_FREE; X_DIMS+1]}, G_, kernel_inputs{:});
            Gnext_(goal_grid) = 0;
            Gaction = cost_fixed + cost_action_ + gamma_*Gnext_;
            Gaction(goal_grid) = 0;
            
            actionUpdate = Gaction < minG;
            actionUpdateBar = ~actionUpdate;
            minG = (Gnext_ .* actionUpdate) + (minG .* actionUpdateBar);
            for uui = 1:1:length(U_DIMS_FREE)
                uu = U_DIMS_FREE(uui);
                unew{uu} = (u_{uu} .* actionUpdate) + (unew{uu} .* actionUpdateBar);
                unew{uu}(goal_grid) = u0(uu);
            end
        end

        uerror_max = gpuArray(zeros(length(U_DIMS_FREE), 1));
        uerror_mean = gpuArray(zeros(length(U_DIMS_FREE), 1));
        for uui = 1:1:length(U_DIMS_FREE)
            uu = U_DIMS_FREE(uui);
            uerror = abs(unew{uu} - u{uu});
            uerror_max(uui) = max(uerror, [], 'all');
            uerror_mean(uui) = sqrt(sum(uerror.^2, 'all')) / numel(G_);
        end

        if (all(uerror_max < u_max_tol(U_DIMS_FREE)) ...
            && all(uerror_mean < u_mean_tol(U_DIMS_FREE)))
            constant_policy_count = constant_policy_count + 1;
            if (constant_policy_count==50)
                fprintf('Policy converged!\n')
                break;
            end
            disp(strcat('Constant Policy count : ', num2str(constant_policy_count)));
        else
            disp(strcat(sprintf('Iteration : %d,', iter), ...
                        ' Max Err :', sprintf(' %.3f', uerror_max), ...
                        ', Mean Err :', sprintf(' %.5f', uerror_mean)));
            constant_policy_count = 0;
        end

        u = unew;
        time_policy_update = time_policy_update + toc(time_policy_update_start);
        clear 'x__' 'u_' 'minG' 'cost_action_' 'Gnext_' 'actionUpdate' 'actionUpdateBar' 'unew';
    end
    time_total = toc(time_start);

%% Return
    policy = u(U_DIMS_FREE);
    info.state_grid = x([X_DIMS_FREE; X_DIMS+1]);
    info.value = G_;
    info.iter = iter;
    info.time_total = time_total;
    info.time_policy_eval = time_policy_eval;
    info.time_policy_update = time_policy_update;

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
