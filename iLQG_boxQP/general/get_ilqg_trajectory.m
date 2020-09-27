function [X, k, K, sub_trajectories_close, ...
          ilqg_cost, ilqg_trace, ilqg_time, ...
          Xinit, Uinit, Costinit] = get_ilqg_trajectory(sys, Op, starts, K_LQR, ...
                                                        sub_policies_LQR, sub_policies_DDP)
%%  Initialize
    
    NUM_CTRL = round(sys.T / sys.dt);
    u0 = sys.u0(sys.U_DIMS_FREE);
    lims = sys.lims(sys.U_DIMS_FREE, :);
    l_point = sys.l_point(sys.X_DIMS_FREE);
    
    Xinit = nan(length(sys.X_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    Uinit = nan(length(sys.U_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    Costinit = zeros(1, size(starts, 2));
    
    if (isfield(sys, 'u0init') && (sys.u0init))
        for kk=1:1:size(starts,2)
            Xinit(:,1,kk) = starts(:,kk);
            discount = 1;
            for ii=1:1:NUM_CTRL
                Uinit(:,ii,kk) = u0;
                Xinit(:,ii+1,kk) = dyn_subs_finite(sys, Xinit(:,ii,kk), Uinit(:,ii,kk), sub_policies_LQR, sys.dt);
                if (any(isnan(Uinit(:,ii,kk))) || any(isnan(Xinit(:,ii,kk))))
                    disp('Check X, U FFull');
                end
                Costinit(kk) = Costinit(kk) + discount*cost_subs(sys, Xinit(:,ii,kk), Uinit(:,ii,kk), sub_policies_LQR)*sys.dt;
                discount = discount*sys.gamma_;
            end
            Costinit(kk) = Costinit(kk) + discount*cost_subs(sys, Xinit(:,NUM_CTRL+1,kk), u0, sub_policies_LQR)*sys.dt;
        end
    else
        for kk=1:1:size(starts,2)
            Xinit(:,1,kk) = starts(:,kk);
            discount = 1;
            for ii=1:1:NUM_CTRL
                Uinit(:,ii,kk) = min(max(u0 + K_LQR*(Xinit(:,ii,kk) - l_point), ...
                                      lims(:,1)), ...
                                  lims(:,2));
                Xinit(:,ii+1,kk) = dyn_subs_finite(sys, Xinit(:,ii,kk), Uinit(:,ii,kk), sub_policies_LQR, sys.dt);
                if (any(isnan(Uinit(:,ii,kk))) || any(isnan(Xinit(:,ii,kk))))
                    disp('Check X, U FFull');
                end
                Costinit(kk) = Costinit(kk) + discount*cost_subs(sys, Xinit(:,ii,kk), Uinit(:,ii,kk), sub_policies_LQR)*sys.dt;
                discount = discount*sys.gamma_;
            end
            Costinit(kk) = Costinit(kk) + discount*cost_subs(sys, Xinit(:,NUM_CTRL+1,kk), u0, sub_policies_LQR)*sys.dt;
        end
    end
    
%%  Compute DDP Trajectory
    
    ilqg_system = @(x, u, sub_policies_, i) ...
                    system_cst(sys, x, u, sub_policies_, sys.full_DDP);
    ilqg_time = zeros(1, size(starts,2));
    ilqg_cost = zeros(1, size(starts,2));
    ilqg_trace = zeros(1, size(starts,2));
    X = zeros(length(sys.X_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    k = zeros(length(sys.U_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    K = zeros(length(sys.U_DIMS_FREE), length(sys.X_DIMS_FREE), NUM_CTRL+1, size(starts, 2));
    sub_trajectories_close = cell(size(sub_policies_DDP));
    sub_trajectories_close(:, 1:2) = sub_policies_DDP(:, 1:2);
    
    for kk=1:1:size(starts, 2)
        
        sub_policies_DDP_ = cell(size(sub_policies_DDP));
        sub_policies_DDP_(:, 1:2) = sub_policies_DDP(:, 1:2);
        sub_policies_DDP_(:, 3) = cellfun(@(x) x(:,:, kk), ...
                                          sub_policies_DDP(:, 3), ...
                                          'UniformOutput', false);
        sub_policies_DDP_(:, 4) = cellfun(@(x) x(:,:,:, kk), ...
                                          sub_policies_DDP(:, 4), ...
                                          'UniformOutput', false);
        sub_policies_DDP_(:, 5) = cellfun(@(x) x(:,:, kk), ...
                                          sub_policies_DDP(:, 5), ...
                                          'UniformOutput', false);
%         Op.cost = Costinit(kk);
        tic;
        [X(:,:,kk), k(:,1:NUM_CTRL,kk), K(:,:,1:NUM_CTRL,kk), ...
         sub_trajectories_close_, ~, ~, ilqg_cost, ilqg_trace] = iLQGGeneral(ilqg_system, ...
                                                                            Xinit(:,1,kk), ... % Don't pass the state trajectory in initialization
                                                                            Uinit(:,1:NUM_CTRL,kk), ...
                                                                            sub_policies_DDP_, ...
                                                                            Op);
        ilqg_time = toc;
        sub_trajectories_close(:,3) = cellfun(@(x, y) cat(3, x, y), ...
                                              sub_trajectories_close(:,3), sub_trajectories_close_(:,3), ...
                                              'UniformOutput', false);
        sub_trajectories_close(:,4) = cellfun(@(x, y) cat(4, x, y), ...
                                              sub_trajectories_close(:,4), sub_trajectories_close_(:,4), ...
                                              'UniformOutput', false);
        sub_trajectories_close(:,5) = cellfun(@(x, y) cat(3, x, y), ...
                                              sub_trajectories_close(:,5), sub_trajectories_close_(:,5), ...
                                              'UniformOutput', false);
        k(:,NUM_CTRL+1,kk) = u0;
        K(:,:,NUM_CTRL+1,kk) = zeros(length(sys.U_DIMS_FREE), length(sys.X_DIMS_FREE));
        disp('Final init point : ');
        Xinit(:, end, kk)
        disp('Final DDP point : ');
        X(:, end, kk)
    end
    
%     save(strcat('XU', sprintf('%d', sys.X_DIMS_FREE), '.mat'), ...
%          'sys', 'X', 'k', 'K', 'u0', 'l_point', 'Xinit', 'Uinit', 'Costinit');
end