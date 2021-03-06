clear;
close all;
clc;

%% 

load('data/Biped2DSystem.mat');

% Debug Mutation Samples
num_samples = 500;
population = generate_population(sys, num_samples);
population_mutated = mutationfunction(sys, linspace(1,num_samples,num_samples), [], [], [], [], [], population);
satisfy_constraint_before = zeros(num_samples,1);
satisfy_constraint_after = zeros(num_samples,1);
for ii=1:1:num_samples
    ii
    p_before = reshape(population(ii, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s_before = reshape(population(ii, (2*sys.U_DIMS+1):end), sys.U_DIMS, sys.X_DIMS);
    [c_before, ceq_before] = constraints(p_before, s_before);
    satisfy_constraint_before(ii) = (all(c_before <= 0) && all(ceq_before==0));
    p_after = reshape(population_mutated(ii, 1:(2*sys.U_DIMS)), sys.U_DIMS, 2);
    s_after = reshape(population_mutated(ii, (2*sys.U_DIMS+1):end), sys.U_DIMS, sys.X_DIMS);
    [c_after, ceq_after] = constraints(p_after, s_after);
    satisfy_constraint_after(ii) = (all(c_after <= 0) && all(ceq_after==0));
end

if (sum(satisfy_constraint_before) ~= num_samples)
    disp('Check generate population');
    return;
end    

if (sum(satisfy_constraint_after) == num_samples)
    disp('All mutations are valid');
else
    invalid_mutations = find(satisfy_constraint_after==0);
    disp('Invalid mutations can be found in invalid_mutations');
end

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