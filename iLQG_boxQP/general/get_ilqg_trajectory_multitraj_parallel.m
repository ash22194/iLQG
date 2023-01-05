function [XFinal, UFinal, X, k, K, sub_trajectories_close, ...
          ilqg_cost, ilqg_trace, ilqg_time, ...
          Xinit, Uinit, Costinit] = get_ilqg_trajectory_multitraj_parallel(sys, Op, starts, K_LQR, ...
                                                        sub_policies_LQR, sub_policies_DDP)
%%  Initialize
    
    NUM_CTRL = round(sys.T / sys.dt);
    u0 = sys.u0(sys.U_DIMS_FREE);
    lims = sys.lims(sys.U_DIMS_FREE, :);
    l_point = sys.l_point(sys.X_DIMS_FREE);
    
    Xinit = nan(length(sys.X_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    Uinit = nan(length(sys.U_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    Costinit = zeros(1, size(starts, 2));
    
    Mdl = cellfun(@(x) KDTreeSearcher(x'), sub_policies_DDP(:,6), 'UniformOutput', false);
    
    if (isfield(sys, 'u0init') && (sys.u0init==1))
        for kk=1:1:size(starts,2)
            [Xinit(:,:,kk), Uinit(:,:,kk), Costinit(kk)] = rollout_trajectory(sys, NUM_CTRL, starts(:,kk), ...
                                               u0, zeros(length(u0), length(sys.X_DIMS_FREE)), sub_policies_DDP, Mdl);
        end
    elseif (isfield(sys, 'u0init') && (sys.u0init==0))
        for kk=1:1:size(starts,2)
            [Xinit(:,:,kk), Uinit(:,:,kk), Costinit(kk)] = rollout_trajectory(sys, NUM_CTRL, starts(:,kk), ...
                                                                                    u0, K_LQR, sub_policies_DDP, Mdl);
        end
    else
        % Pick the least cost initial trajectory
        parfor kk=1:1:size(starts,2)
            [Xinit0, Uinit0, Costinit0] = rollout_trajectory(sys, NUM_CTRL, starts(:,kk), ...
                                               u0, zeros(length(u0), length(sys.X_DIMS_FREE)), sub_policies_DDP, Mdl);
            [XinitLQR, UinitLQR, CostinitLQR] = rollout_trajectory(sys, NUM_CTRL, starts(:,kk), ...
                                               u0, K_LQR, sub_policies_DDP, Mdl);
            if (Costinit0 < CostinitLQR)
                Xinit(:,:,kk) = Xinit0;
                Uinit(:,:,kk) = Uinit0;
                Costinit(kk) = Costinit0;
            else
                Xinit(:,:,kk) = XinitLQR;
                Uinit(:,:,kk) = UinitLQR;
                Costinit(kk) = CostinitLQR;
            end
        end
    end

%%  Compute DDP Trajectory
    
    ilqg_system = @(x, u, sub_policies_, i) ...
                    system_cst(sys, x, u, sub_policies_, sys.full_DDP);
    ilqg_time = zeros(1, size(starts,2));
    ilqg_cost = zeros(1, NUM_CTRL+1, size(starts,2));
    ilqg_trace = cell(size(starts,2), 1);
    X = zeros(length(sys.X_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    XFinal = zeros(length(sys.X_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    k = zeros(length(sys.U_DIMS_FREE), NUM_CTRL, size(starts, 2));
    UFinal = zeros(length(sys.U_DIMS_FREE), NUM_CTRL, size(starts, 2));
    K = zeros(length(sys.U_DIMS_FREE), length(sys.X_DIMS_FREE), NUM_CTRL, size(starts, 2));
    sub_trajectories_close = cell(size(sub_policies_DDP));
    sub_trajectories_close(:, 1:2) = sub_policies_DDP(:, 1:2);
    sub_trajectories_close_ = cell(size(starts, 2), 1);
    
    parfor kk=1:1:size(starts, 2)
        disp(strcat(num2str(kk),') Final init point : ', sprintf('%.4f ', Xinit(:, end, kk))));
%         Op.cost = Costinit(kk);
        tic;
        [XFinal_, X(:,:,kk), UFinal(:,:,kk), k(:,:,kk), K(:,:,:,kk), ...
         sub_trajectories_close_{kk, 1}, ~, ~, ilqg_cost(:,:,kk), ilqg_trace{kk}] = iLQGGeneralKDTree(ilqg_system, ...
                                                                            Xinit(:,1,kk), ... % Don't pass the state trajectory in initialization
                                                                            Uinit(:,1:NUM_CTRL,kk), ...
                                                                            sub_policies_DDP, ...
                                                                            Op);
        ilqg_time(kk) = toc;
        XFinal(:,:,kk) = XFinal_;
        disp(strcat(num2str(kk),') Final init point : ', sprintf('%.4f ', Xinit(:, end, kk))));
        disp(strcat(num2str(kk),') Final DDP point : ', sprintf('%.4f ', XFinal_(:, end)), sprintf('\n')));
    end
    
    k = cat(2, k, repmat(u0, [1, 1, size(starts, 2)]));
    UFinal = cat(2, UFinal, repmat(u0, [1, 1, size(starts, 2)]));
    K = cat(3, K, zeros(length(sys.U_DIMS_FREE), length(sys.X_DIMS_FREE), 1, size(starts, 2)));
    for kk=1:1:size(starts, 2)
        
        sub_trajectories_close(:,3) = cellfun(@(x, y) cat(3, x, y), ...
                                              sub_trajectories_close(:,3), sub_trajectories_close_{kk,1}(:,3), ...
                                              'UniformOutput', false);
        sub_trajectories_close(:,4) = cellfun(@(x, y) cat(4, x, y), ...
                                              sub_trajectories_close(:,4), sub_trajectories_close_{kk,1}(:,4), ...
                                              'UniformOutput', false);
        sub_trajectories_close(:,5) = cellfun(@(x, y) cat(3, x, y), ...
                                              sub_trajectories_close(:,5), sub_trajectories_close_{kk,1}(:,5), ...
                                              'UniformOutput', false);
        sub_trajectories_close(:,6) = cellfun(@(x, y) cat(3, x, y), ...
                                              sub_trajectories_close(:,6), sub_trajectories_close_{kk,1}(:,6), ...
                                              'UniformOutput', false);
    end
    
    X = reshape(X, size(X,1), size(X,2)*size(X,3));
    XFinal = reshape(XFinal, size(XFinal,1), size(XFinal,2)*size(XFinal,3));
    k = reshape(k, size(k,1), size(k,2)*size(k,3));
    UFinal = reshape(UFinal, size(UFinal,1), size(UFinal,2)*size(UFinal,3));
    K = reshape(K, size(K,1), size(K,2), size(K,3)*size(K,4));
    
end

function [X, U, C] = rollout_trajectory(sys, NUM_CTRL, start, u0, K0, sub_policies, Mdl)
    X = nan(length(sys.X_DIMS_FREE), NUM_CTRL+1);
    U = nan(length(sys.U_DIMS_FREE), NUM_CTRL+1);
    C = 0;
    sub_policies_ = cell(size(sub_policies));
    sub_policies_(:,1:2) = sub_policies(:,1:2);

    X(:, 1) = start;
    discount = 1;
    for ii=1:1:NUM_CTRL
        U(:,ii) = u0 + K0*(X(:,ii) - sys.l_point(sys.X_DIMS_FREE));

        for jj=1:1:size(sub_policies, 1)
            closest_x = knnsearch(Mdl{jj}, X(sub_policies{jj,2}, ii)');
            sub_policies_{jj, 3} = sub_policies{jj, 3}(:, closest_x);
            sub_policies_{jj, 4} = sub_policies{jj, 4}(:,:, closest_x);
            sub_policies_{jj, 5} = sub_policies{jj, 5}(:, closest_x);
        end

        X(:,ii+1) = dyn_subs_finite2(sys, X(:,ii), U(:,ii), sub_policies_, sys.dt);
        if (any(isnan(U(:,ii))) || any(isnan(X(:,ii))))
            disp('Check X, U FFull');
        end
        
        C = C + discount*cost_subs(sys, X(:,ii), U(:,ii), sub_policies_)*sys.dt;
        discount = discount*sys.gamma_;
    end
    for jj=1:1:size(sub_policies, 1)
        closest_x = knnsearch(Mdl{jj}, X(sub_policies{jj,2}, NUM_CTRL+1)');
        sub_policies_{jj, 3} = sub_policies{jj, 3}(:, closest_x);
        sub_policies_{jj, 4} = sub_policies{jj, 4}(:,:, closest_x);
        sub_policies_{jj, 5} = sub_policies{jj, 5}(:, closest_x);
    end
    C = C + discount*cost_subs(sys, X(:,NUM_CTRL+1), u0, sub_policies_)*sys.dt;
end