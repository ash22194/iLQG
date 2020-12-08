clear;
close all;
clc;

%% 

restoredefaultpath;
n = 3;
system_name = sprintf('manipulator%ddof', n);
addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
addpath(genpath('test'));
load(strcat('../iLQG_boxQP/new_systems/', system_name, '/sys.mat'), 'sys');

sys.X_DIMS = 2*sys.n; % [thi, ... dthi, ...]
sys.U_DIMS = sys.n;   % [taui]
if (n==2)
    sys.m = [2.5; 0.5]/2; % kg
    sys.l = [0.5; 0.25]/2; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8, 8, 0.6, 0.6])/5;
    sys.R = diag(0.003*(Izz(1)./Izz).^2);
    sys.lims = 5*[-Izz/Izz(1), Izz/Izz(1)]; % action limits
    sys.num_points = 31 * ones(1, sys.X_DIMS);
    sys.num_action_samples = [15, 5];
    
    % Define decompositions to test
    u_x = [];
    pseudo_inputs{1,1} = {{1, 2}}; % r = 2
    
    state_assignments{1} = {{1, 2}; {[], [1,2]}};
    
elseif (n==3)
    sys.m = [2.5; 0.5; 0.1] * 1.1; % kg
    sys.l = [0.5; 0.25; 0.125]; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8*ones(1,3), 0.6*ones(1,3)])/5;
    sys.R = diag(0.004*(Izz(1)./Izz));
    sys.lims = [-16, 16; -7.5, 7.5; -1, 1]; % action limits
    sys.num_points = 31 * ones(1, sys.X_DIMS);
    sys.num_action_samples = [15, 5, 2];
    
    % Define decompositions to test
    u_x = [];
    pseudo_inputs{1,1} = {{1, [2,3]}; {2, [1,3]}; {3, [1,2]}}; % r = 2
    pseudo_inputs{2,1} = {{1, 2, 3}}; % r = 3
    
    state_assignments{1} = {{1, [2,3]}; {2, [1,3]}; {3, [1,2]}; {[], [1,2,3]}};
    state_assignments{2} = {{1, 2, 3}; {1, [2,3], []}; {2, [1,3], []}; {3, [1,2], []}; {[], [1,2,3], []}};

elseif (n==4)
    sys.m = [5.4; 1.8; 0.6; 0.2]; % kg
    sys.l = [0.2; 0.5; 0.25; 0.125]; % m
    Izz = sys.m.*((sys.l));
    sys.Q = diag([8*ones(1,4), 0.2*ones(1,4)])/2;
    sys.R = diag([0.002; 0.004*(Izz(2)./Izz(2:end))]);
    sys.lims = [-24, 24; -15, 15; -7.5, 7.5; -1, 1]; % action limits
    sys.num_points = 31 * ones(1, sys.X_DIMS);
    sys.num_action_samples = [24, 15, 8, 2];
    
    % Define decompositions to test
    u_x = [];
    pseudo_inputs{1,1} = {{1, [2,3,4]}; {2, [1,3,4]}; {3, [1,2,4]}; {4, [1,2,3]};
                        {[1,2], [3,4]}; {[1,3], [2,4]}; {[1,4], [2,3]}}; % r = 2
    pseudo_inputs{2,1} = {{1, 2, [3,4]}; {1, 3, [2,4]}; {1, 4, [2,3]}; {2, 3, [1,4]}; {2, 4, [1,3]}; {3, 4, [1,2]}}; % r = 3
    pseudo_inputs{3,1} = {{1, 2, 3, 4}}; % r = 4
    
    state_assignments{1} = {{1, [2,3,4]}; {2, [1,3,4]}; {3, [1,2,4]}; {4, [1,2,3]};
                            {[1,2], [3,4]}; {[1,3], [2,4]}; {[1,4], [2,3]}; {[], [1,2,3,4]}};
    state_assignments{2} = {{1, [2,3,4], []}; {2, [1,3,4], []}; {3, [1,2,4], []}; {4, [1,2,3], []};
                            {[1,2], [3,4], []}; {[1,3], [2,4], []}; {[1,4], [2,3], []}; {[], [1,2,3,4], []};
                            {1, 2, [3,4]}; {1, 3, [2,4]}; {1, 4, [2,3]}; {2, 3, [1,4]}; {2, 4, [1,3]}; {3, 4, [1,2]}};
    state_assignments{3} = {{1, [2,3,4], [], []}; {2, [1,3,4], [], []}; {3, [1,2,4], [], []}; {4, [1,2,3], [], []};
                            {[1,2], [3,4], [], []}; {[1,3], [2,4], [], []}; {[1,4], [2,3], [], []}; {[], [1,2,3,4], [], []};
                            {1, 2, [3,4], []}; {1, 3, [2,4], []}; {1, 4, [2,3], []}; {2, 3, [1,4], []}; {2, 4, [1,3], []}; {3, 4, [1,2], []};
                            {1, 2, 3, 4}};
end

for r=1:1:size(pseudo_inputs, 1)
    for sub_r=1:1:size(pseudo_inputs{r}, 1)
        current_pseudo_inputs = pseudo_inputs{r}{sub_r};
        orders_of_inputs = perms(1:1:(r+1));
        for ori=1:1:size(orders_of_inputs,1)
            input_order = orders_of_inputs(ori,:);
            ordered_pseudo_input = current_pseudo_inputs(input_order);
            % Define action tree for purely cascaded
            p_casc = [zeros(n,1), ones(n,1)];
            for ori_=2:1:(r+1)
                p_casc(ordered_pseudo_input{ori_}, 1) = ordered_pseudo_input{ori_-1}(1);
            end

            for sub_s=1:1:size(state_assignments{r}, 1)
                current_state_assignment = state_assignments{r}{sub_s};
                orders_of_states = perms(1:1:(r+1));
                perms_of_states = [];

                for ors=1:1:size(orders_of_states, 1)
                   state_order = orders_of_states(ors,:);
                   ordered_states = current_state_assignment(state_order);
                   zero_ordered_states = cell2mat(cellfun(@(x) zero_empty(x), ordered_states, 'UniformOutput', false));
                   perm_id = sum(zero_ordered_states.*((n+1).^(linspace(0, length(zero_ordered_states)-1, length(zero_ordered_states)))));
                   if (any(perms_of_states == perm_id))
                       continue;
                   else
                       perms_of_states = [perms_of_states; perm_id];
                   end
                   if (isempty(ordered_states{r+1}))
                       continue;
                   end
                   % Define state assignment
                   s = zeros(n, 2*n);
                   for ors_=1:1:(r+1)
                       s(ordered_pseudo_input{ors_}, ordered_states{ors_}) = 1;
                       s(ordered_pseudo_input{ors_}, ordered_states{ors_} + n) = 1;
                   end
                   u_x = [u_x; reshape(p_casc, 1, 2*n), reshape(s, 1, 2*n^2)];
                end
            end
        end

        % Define action tree for purely decoupled
        p_dec = zeros(n,2);
        for ori_=1:1:(r+1)
            p_dec(current_pseudo_inputs{ori_}, 2) = ori_;
        end

        for sub_s=1:1:size(state_assignments{r}, 1)
            current_state_assignment = state_assignments{r}{sub_s};
            orders_of_states = perms(1:1:(r+1));
            perms_of_states = [];

            for ors=1:1:size(orders_of_states, 1)
               state_order = orders_of_states(ors,:);
               ordered_states = current_state_assignment(state_order);
               zero_ordered_states = cell2mat(cellfun(@(x) zero_empty(x), ordered_states, 'UniformOutput', false));
               perm_id = sum(zero_ordered_states.*((n+1).^(linspace(0, length(zero_ordered_states)-1, length(zero_ordered_states)))));
               if (any(perms_of_states == perm_id))
                   continue;
               else
                   perms_of_states = [perms_of_states; perm_id];
               end
               if any(cellfun(@(x) isempty(x), ordered_states))
                   continue;
               end
               % Define state assignment
               s = zeros(n, 2*n);
               for ors_=1:1:(r+1)
                   s(ordered_pseudo_input{ors_}, ordered_states{ors_}) = 1;
                   s(ordered_pseudo_input{ors_}, ordered_states{ors_} + n) = 1;
               end
               u_x = [u_x; reshape(p_dec, 1, 2*n), reshape(s, 1, 2*n^2)];
            end
        end
    end
end

sys.max_iter = 2000;
sys.max_policy_iter = 100;
sys.g = 9.81; % m/s^2
sys.dt = 0.001;
sys.gamma_ = 0.997;
sys.lambda_ = (1 - sys.gamma_) / sys.dt;

sys.l_point = zeros(sys.X_DIMS, 1);
sys.l_point(1) = pi;
sys.goal = sys.l_point;
sys.u0 = zeros(sys.U_DIMS, 1);
sys.fxfu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];

fxfu = sys.fxfu_func(sys.l_point, sys.u0);
sys.A = fxfu(:,1:sys.X_DIMS);
sys.B = fxfu(:,(1+sys.X_DIMS):end);
[~, S_joint, ~] = lqr(sys.A - eye(size(sys.A,1))*sys.lambda_/2, sys.B, sys.Q, sys.R, zeros(size(sys.A,1), size(sys.B,2)));
sys.S =  sym('S', [sys.X_DIMS, sys.X_DIMS]);
sys.S = tril(sys.S,0) + tril(sys.S,-1).';
sys.a = sym('a', [sys.X_DIMS, 1]);
sys.b = sym('b', [sys.X_DIMS, 1]);
sys.err_lqr = sys.x.'*sys.S*sys.x - sys.x.'*S_joint*sys.x;
for ii=1:1:sys.X_DIMS
    sys.err_lqr = int(sys.err_lqr, sys.x(ii), [sys.a(ii), sys.b(ii)]);
end
sys.err_lqr_func = matlabFunction(simplify(sys.err_lqr), 'Vars', {sys.S, sys.a, sys.b});

sys.state_bounds = [repmat([-pi/3, pi/3], [n, 1]);
                    repmat([-0.5, 0.5], [n, 1])];
sys.state_bounds(1,:) = sys.state_bounds(1,:) + pi;
sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1));

err_lqr = zeros(1, size(u_x,1));
compute_fraction = zeros(1, size(u_x,1));
for dd=1:1:size(u_x,1)
    p = reshape(u_x(dd, 1:(2*sys.U_DIMS)), [sys.U_DIMS, 2]);
    s = reshape(u_x(dd, (1+2*sys.U_DIMS):end), [sys.U_DIMS, sys.X_DIMS]);
    err_lqr(dd) = computeLQRMeasure(sys, p, s);
    compute_fraction(dd) = computeComplexityEstimates(sys, p, s);
end

pareto_front = paretoQS([err_lqr', compute_fraction']);
pareto_front(~(err_lqr(pareto_front) < inf)) = []; % Remove candidates with infinite LQR error
u_xp = u_x(pareto_front, :);
% save(strcat('data/', system_name, '_paretofront.mat'), 'u_xp');

%% Functions

function x = zero_empty(x)
    if (isempty(x))
       x = 0;
    end
end
