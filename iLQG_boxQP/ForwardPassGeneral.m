function [X, U, c] = ForwardPassGeneral(sys, sub_policies, x_start)

    NUM_STARTS = size(x_start, 2);
    NUM_CTRL = round(sys.T / sys.dt);
    X_DIMS = sys.X_DIMS;
    U_DIMS = sys.U_DIMS;

    X = zeros(X_DIMS, NUM_CTRL + 1, NUM_STARTS);
    U = zeros(U_DIMS, NUM_CTRL, NUM_STARTS);
    c = zeros(1, NUM_STARTS);
    
    for kk = 1:1:NUM_STARTS
        
        fprintf('Trajectory : %d / %d\n', kk, NUM_STARTS);
        
        X(:, 1, kk) = x_start(:, kk);
        discount = 1;
        for tt = 1:1:NUM_CTRL
            closest_x = cellfun(@(x, y) min_index(vecnorm(X(x, tt, kk) - y(:,:, kk), 2, 1)), sub_policies(:, 2), sub_policies(:, 5));
            for jj=1:1:size(sub_policies, 1)
                U(sub_policies{jj, 1}, tt, kk) = sub_policies{jj, 3}(:, closest_x(jj), kk) ...
                                                 + sub_policies{jj, 4}(:,:, closest_x(jj), kk) ...
                                                   * (X(sub_policies{jj, 2}, tt, kk) - sub_policies{jj, 5}(:, closest_x(jj), kk));
            end
            U(:, tt, kk) = max(sys.lims(:,1), min(sys.lims(:,2), U(:, tt, kk)));
            X(:, tt+1, kk) = dyn_finite(sys, X(:, tt, kk), U(:, tt, kk), sys.dt);
            c(kk) = c(kk) + discount * cost(sys, X(:, tt, kk), U(:, tt, kk)) * sys.dt;
            discount = discount * sys.gamma_;
        end
    end

end

function out = min_index(x)
    [~, out] = min(x);
end