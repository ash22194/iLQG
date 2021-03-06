clear;
close all;
clc;

%%
load('data/QuadcopterSystem.mat');
sys.xu = [sys.x; sys.u];
d = clock;
rng(d(end));

num_samples = 10600;
population = generate_population(sys, num_samples);
err_lqr = zeros(num_samples, 1);
for ii=1:num_samples
    ii
    p = reshape(population(ii, 1:2*sys.U_DIMS), sys.U_DIMS, 2);
    s = reshape(population(ii, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);
    err_lqr(ii) = computeLQRMeasure(sys, p, s);
end

[~, unique_population_id] = extract_decompositions_from_population(sys, population);
unique_population = population(unique_population_id,:);
err_lqr = err_lqr(unique_population_id);

[err_lqr_sorted, err_lqr_sorted_id] = sort(err_lqr);

[best_score, best_decomposition] = min(err_lqr);
p = reshape(unique_population(best_decomposition, 1:2*sys.U_DIMS), sys.U_DIMS, 2);
s = reshape(unique_population(best_decomposition, (1+2*sys.U_DIMS):end), sys.U_DIMS, sys.X_DIMS);

%% Functions

function population = generate_population(sys, n)
    
    % Decide input coupling
    invalid = true(1, n);
    while (sum(invalid) > 0)
       r(:, invalid) = randi([1,sys.U_DIMS], sys.U_DIMS, sum(invalid));
       invalid = vecnorm(r - mean(r, 1)) < 1e-4;
    end
    p = zeros(sys.U_DIMS, 2, n);
    s = zeros(sys.U_DIMS, sys.X_DIMS, n);
    
    for ii=1:1:n
        rC = unique(r(:,ii));
        pseudo_inputs = cell(length(rC)+1, 1);
        pseudo_inputs{end, 1} = [0];
        
        % Random state assignments for the inputs
        state = linspace(1, sys.X_DIMS, sys.X_DIMS);
        for jj=1:1:length(rC)
            pseudo_inputs{jj} = find(r(:,ii)==rC(jj));
            
            k = randi([1, length(state)-(length(rC) - jj)], 1);
            sub_state = nchoosek(linspace(1,length(state),length(state)), k);
            sub_state = sub_state(randi([1, size(sub_state,1)],1),:);
            sub_state_ = state(sub_state);
            state(sub_state) = [];

            s(r(:,ii)==rC(jj), sub_state_, ii) = 1;
        end
        s(r(:,ii)==rC(jj), state, ii) = 1;
        
        % Use random Prufer codes to generate a tree
        prufer_seq = randi([1, length(rC) + 1], length(rC)-1, 1);
        node_list = linspace(1, length(rC) + 1, length(rC) + 1)';
        child_count = zeros(length(rC) + 1, 1);
        while(length(node_list) > 2)
            child = min(setdiff(node_list, prufer_seq));
            parent = prufer_seq(1);
            child_count(parent) = child_count(parent) + 1;
            
            node_list(node_list == child) = [];
            prufer_seq = prufer_seq(2:end);
            
            p(pseudo_inputs{child}, 1, ii) = pseudo_inputs{parent}(1);
            p(pseudo_inputs{child}, 2, ii) = child_count(parent);
        end
        child = node_list(1);
        parent = node_list(2);
        child_count(parent) = child_count(parent) + 1;
        p(pseudo_inputs{child}, 1, ii) = pseudo_inputs{parent}(1);
        p(pseudo_inputs{child}, 2, ii) = child_count(parent);
        
    end
    
    population = [reshape(p, 2*sys.U_DIMS, n)', reshape(s, sys.U_DIMS*sys.X_DIMS, n)'];

end