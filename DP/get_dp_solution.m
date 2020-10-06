function [policy, value, iter] = get_dp_solution(sys, Op, sub_policies)
%% Parameters
    X_DIMS            = sys.X_DIMS;
    X_DIMS_FREE       = sys.X_DIMS_FREE;
    X_DIMS_FIXED      = sys.X_DIMS_FIXED;
    U_DIMS            = sys.U_DIMS;
    U_DIMS_FREE       = sys.U_DIMS_FREE;
    U_DIMS_FIXED      = sys.U_DIMS_FIXED;
    U_DIMS_CONTROLLED = sys.U_DIMS_CONTROLLED;
    limits            = sys.limits; % state limits
    lims              = sys.lims;   % action limits
    goal              = sys.goal;
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
    num_action_samples = prod(Op.num_action_samples(U_DIMS_FREE));
    
%% Initialize 
    % Create state grids
    x = cell(X_DIMS, 1);
    x(X_DIMS_FIXED) = num2cell(goal(X_DIMS_FIXED));
    grid_indices = cell(length(X_DIMS_FREE), 1);
    for xxi = 1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        grid_indices{xxi} = linspace(limits(xx,1), limits(xx,2), num_points(xx));
    end
    [x{X_DIMS_FREE}] = ndgrid(grid_indices{:});
    
    % Initialize policy grids - fixed, free and controlled
    u = cell(U_DIMS, 1);
    u(U_DIMS_FIXED) = num2cell(zeros(length(U_DIMS_FIXED), 1));
    for uui = 1:1:length(U_DIMS_FREE)
        uu = U_DIMS_FREE(uui);
        u{uu} = lims(uu, 1) + (lims(uu, 2) - lims(uu, 1)) * rand(num_points(X_DIMS_FREE));
    end
    for jj=1:1:size(sub_policies, 1)
        U_SUBDIM = sub_policies{jj,1};
        X_SUBDIM = sub_policies{jj,2};
        X_SUBDIM_BAR = 1:X_DIMS;
        X_SUBDIM_BAR(X_SUBDIM) = [];
        
        subpolicy_size = num_points;
        subpolicy_size(X_SUBDIM_BAR) = 1;
        subpolicy_newsize = num_points;
        subpolicy_newsize(X_SUBDIM) = 1;
        
        u(U_SUBDIM) = cellfun(@(x) repmat(reshape(x, subpolicy_size), subpolicy_newsize), sub_policies{jj,3}, 'UniformOutput', false);
    end
    
    % Initialize cost functions
    % State cost assuming diagonal Q matrix
    cost_state = 0;
    for xxi = 1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        cost_state = cost_state + dt * Q(xx,xx) * (x{xx} - goal(xx)).^2;
    end
    
    % Controlled action cost assuming diagonal R matrix
    cost_action_controlled = 0;
    for uui=1:1:length(U_DIMS_CONTROLLED)
        uu = U_DIMS_CONTROLLED(uui);
        cost_action_controlled = cost_action_controlled + dt * R(uu, uu) * (u{uu} - u0(uu)).^2;
    end
    cost_fixed = cost_state + cost_action_controlled;
    
    goal_grid = zeros(length(X_DIMS_FREE), 1);
    for xxi=1:1:length(X_DIMS_FREE)
        xx = X_DIMS_FREE(xxi);
        [~, goal_grid(xxi)] = min(abs(grid_indices{xxi} - goal(xx)));
    end
    goal_grid = num2cell(goal_grid);
    
    G_ = zeros(num_points(X_DIMS_FREE));
    F = griddedInterpolant(x{X_DIMS_FREE}, G_);
    constant_policy_count = 0;
    iter = 0;
    
%% Policy Iteration

    while (iter < max_iter)
        iter = iter + 1;
        G = ones(size(G_));
        policy_iter = 0;
        cost_action = 0;
        x_ = dyn_finite_rk4(sys, x, u, dt);
        for uui = 1:1:length(U_DIMS_FREE)
            uu = U_DIMS_FREE(uui);
            cost_action = cost_action + dt * R(uu, uu) * (u{uu} - u0(uu)).^2;
        end
        
        cost_total = cost_fixed + cost_action;
        % Compute value function for current policy
         while ((max(abs(G_ - G), [], 'all') > gtol) && (policy_iter < max_policy_iter))
            % Iterate to estimate value function
            G = G_;
            policy_iter = policy_iter + 1;
            Gnext = F(x_{X_DIMS_FREE});
            Gnext(goal_grid{:}) = 0;
            
            % Update value function
            G_ = cost_total + gamma_*Gnext;
            G_(goal_grid{:}) = 0;
            F.Values = G_;
        end
        clear 'G';

        % Update policy
        % Compute Q(x,u) values g(x,u) + G(x+udt) for different actions when value
        % function G has converged and greedily update the policy i.e. choose
        % action u* at state x such that u* = argmin Q(x,u)
        minG = G_;
        unew = u;
        u_ = u;
        for i = 1:1:num_action_samples
            cost_action_ = 0;
            for uui=1:1:length(U_DIMS_FREE)
                uu = U_DIMS_FREE(uui);
                u_{uu} = lims(uu, 1) + (lims(uu, 2) - lims(uu, 1)) * rand(size(G_));
                cost_action_ = cost_action_ + dt * R(uu, uu) * (u_{uu} - u0(uu)).^2;
            end
            
            x__ = dyn_finite_rk4(sys, x, u_, dt);
            
            Gnext_ = F(x__{X_DIMS_FREE});
            Gnext_(goal_grid{:}) = 0;
            Gaction = cost_fixed + cost_action_ + gamma_*Gnext_;
            Gaction(goal_grid{:}) = 0;
            
            actionUpdate = Gaction < minG; 
            minG(actionUpdate) = Gaction(actionUpdate);
            for uui = 1:1:length(U_DIMS_FREE)
                uu = U_DIMS_FREE(uui);
                unew{uu}(actionUpdate) = u_{uu}(actionUpdate);
                unew{uu}(goal_grid{:}) = u0(uu);
            end
        end
        
        uerror_max = zeros(length(U_DIMS_FREE));
        uerror_mean = zeros(length(U_DIMS_FREE));
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
                        ' Max Err :', sprintf(' %.5f', uerror_max), ...
                        ', Mean Err :', sprintf(' %.5f', uerror_mean)));
            constant_policy_count = 0;
        end
        
        u = unew;
    end
    
%% Return
    policy = u(U_DIMS_FREE);
    value = G_;
end