clear;
close all;
clc;

%% 

system_name = 'manipulator3dof';
addpath('utils');
load(strcat('/usr0/home/akhadke/PolicyDecomposition/MCTS/data/', system_name, 'System.mat'));

sys.X_DIMS_MASKED = eye(sys.X_DIMS);
num_decompositions = calcAllDecompositions(sys.X_DIMS, sys.U_DIMS);

population_new = uniform_input_tree_sampling(sys, 5*num_decompositions);
[~, uid_new] = extract_decompositions_from_population(sys, population_new);

population_count_new = zeros(length(uid_new), 1);
for uui=1:1:length(uid_new)
    population_count_new(uui) = sum(all(population_new(uid_new(uui), :) == population_new, 2));
end
mean_count_new = mean(population_count_new);
std_count_new = std(population_count_new);

population_old = generate_population(sys, 5*num_decompositions);
[~, uid_old] = extract_decompositions_from_population(sys, population_old);

population_count_old = zeros(length(uid_old), 1);
for uui=1:1:length(uid_old)
    population_count_old(uui) = sum(all(population_old(uid_old(uui), :) == population_old, 2));
end
mean_count_old = mean(population_count_old);
std_count_old = std(population_count_old);

%% Functions

function population = generate_population(sys, n)
    
    % Decide input coupling
    invalid = true(1,n);
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
        
        if (any(p(:,1,ii) == linspace(1, sys.U_DIMS, sys.U_DIMS)'))
            disp('Loops!');
        end
        [c, c_eq] = constraints(p(:,:,ii), s(:,:,ii));
        if (any(c_eq~=0) || any(c > 0))
            disp('Invalid sample');
        end
    end
    
    population = [reshape(p, 2*sys.U_DIMS, n)', reshape(s, sys.U_DIMS*sys.X_DIMS, n)'];

end