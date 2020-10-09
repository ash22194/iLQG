function [X, U, c] = ForwardPassGeneral_multitraj(sys, sub_policies, x_start)

    NUM_STARTS = size(x_start, 2);
    NUM_CTRL = round(sys.T / sys.dt);
    X_DIMS = sys.X_DIMS;
    U_DIMS = sys.U_DIMS;

    X = zeros(X_DIMS, NUM_CTRL + 1, NUM_STARTS);
    U = zeros(U_DIMS, NUM_CTRL, NUM_STARTS);
    c = zeros(1, NUM_STARTS);
    
    Mdl = cellfun(@(x) KDTreeSearcher(x'), sub_policies(:,5), 'UniformOutput', false);
    
    for kk = 1:1:NUM_STARTS
        
        fprintf('Trajectory : %d / %d\n', kk, NUM_STARTS);
        
        X(:, 1, kk) = x_start(:, kk);
        discount = 1;
        sub_policies_tt = cell(size(sub_policies));
        sub_policies_tt(:,1:2) = sub_policies(:,1:2);
        for tt = 1:1:NUM_CTRL
%             closest_x = cellfun(@(x, y) min_index(vecnorm(X(x, tt, kk) - y, 2, 1)), sub_policies(:, 2), sub_policies(:, 5));
            for jj=1:1:size(sub_policies, 1)
                closest_x = knnsearch(Mdl{jj}, X(sub_policies_tt{jj, 2}, tt, kk)');
                sub_policies_tt{jj, 3} = sub_policies{jj, 3}(:, closest_x);
                sub_policies_tt{jj, 4} = sub_policies{jj, 4}(:,:, closest_x);
                sub_policies_tt{jj, 5} = sub_policies{jj, 5}(:, closest_x);
            end
            
            X(:,tt+1, kk) = dyn_subs_finite2(sys, X(:, tt, kk), zeros(0,1), sub_policies_tt, sys.dt); 
            c(kk) = c(kk) + discount * cost(sys, X(:, tt, kk), U(:, tt, kk)) * sys.dt;
            discount = discount * sys.gamma_;
        end
    end

end

function out = min_index(x)
    [~, out] = min(x, [], 2);
end