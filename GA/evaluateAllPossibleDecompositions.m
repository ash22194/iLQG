clear;
clc;
close all;

%% 

restoredefaultpath();
system_name = 'biped2d';
addpath(strcat('../iLQG_boxQP/new_systems/', system_name));
addpath('../iLQG_boxQP/iLQG utilities/decomposition_count');
addpath('utils');
load(strcat('../../MCTS/data/', system_name, 'System.mat'), 'sys');

sys.X_DIMS_MASKED = eye(sys.X_DIMS);
% Add function to compute linearized dynamics
if (~isfield(sys, 'fxfu_func'))
    sys.fxfu_func = @(x, u) [dynx(sys, x, u), dynu(sys, x, u)];
    fxfu = sys.fxfu_func(sys.l_point, sys.u0);
    sys.A = fxfu(:, 1:sys.X_DIMS);
    sys.B = fxfu(:, (1+sys.X_DIMS):end);
end

% Add function to compute LQR metric
if (~isfield(sys, 'err_lqr_func'))
    sys.lambda_ = (1 - sys.gamma_) / sys.dt;
    [~, S_j, ~] = lqr(sys.A - (sys.lambda_ / 2) * eye(size(sys.A,1)), sys.B, sys.Q, sys.R, zeros(size(sys.B)));
    sys.S = sym('S', [sys.X_DIMS, sys.X_DIMS]);
    sys.S = tril(sys.S,0) + tril(sys.S,-1).';
    sys.a = sym('a', [sys.X_DIMS, 1]);
    sys.b = sym('b', [sys.X_DIMS, 1]);
    err_lqr = sys.x.'*sys.S*sys.x - sys.x.'*S_j*sys.x;
    for ii=1:1:sys.X_DIMS
        err_lqr = int(err_lqr, sys.x(ii), [sys.a(ii), sys.b(ii)]);
    end
    sys.err_lqr_func = matlabFunction(err_lqr, 'Vars', {sys.S, sys.a, sys.b});
    sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1), 1);
end
sys.measure_func = @(err_lqr, err_compute) (1 - exp(-err_lqr)) .* err_compute;

%% Enumerate Decompositions
number_of_decompositions = calcAllDecompositions(size(sys.X_DIMS_MASKED, 1), sys.U_DIMS);
policy_decompositions = cell(number_of_decompositions, 5);

decomposition_count = 0;

for rr=2:1:sys.U_DIMS
	% Number of nodes in input-tree

	% Possible rr groupings on inputs
	rr_groupings = partitions(sys.U_DIMS, rr);
	assert(size(rr_groupings{1}, 2)==rr);

	for rrg=1:1:size(rr_groupings, 1)
        rr_groupings{rrg} = cat(2, rr_groupings{rrg}, [0]);
		prufer_codes = (rr+1) * ones((rr+1)^(rr-1), rr-1); % sequenc of length (rr - 1) with entries from 1 to (rr + 1)
		for pp = 1:1:size(prufer_codes, 1)
			pp_ = pp;
			for seq = 1:1:(rr-1)
				reminder = mod(pp_, rr+1);
				pp_ = (pp_ - reminder) / (rr+1);
				if (reminder >=1)
					prufer_codes(pp, seq) = reminder;
				end
			end

			% Each prufer code represents an input-tree
			p = zeros(sys.U_DIMS, 2);
			node_list = linspace(1, rr + 1, rr + 1)';
			child_count = zeros(rr + 1, 1);
	        prufer_seq = prufer_codes(pp, :);
	        while(length(node_list) > 2)
	            child = min(setdiff(node_list, prufer_seq));
	            parent = prufer_seq(1);
	            child_count(parent) = child_count(parent) + 1;
	            
	            node_list(node_list == child) = [];
	            prufer_seq = prufer_seq(2:end);
	            
	            p(rr_groupings{rrg}{child}, 1) = rr_groupings{rrg}{parent}(1);
	            p(rr_groupings{rrg}{child}, 2) = child_count(parent);
	        end
	        child = node_list(1);
	        parent = node_list(2);
	        child_count(parent) = child_count(parent) + 1;
	        p(rr_groupings{rrg}{child}, 1) = rr_groupings{rrg}{parent}(1);
	        p(rr_groupings{rrg}{child}, 2) = child_count(parent);

	        % State assignments
	        node_list = linspace(1, rr + 1, rr + 1)';
	        leaf_node_list = linspace(1, rr + 1, rr + 1)';
	        leaf_node_list([prufer_codes(pp, :), (rr + 1)]) = [];
	        non_leaf_node_list = linspace(1, rr + 1, rr + 1)';
	        non_leaf_node_list(leaf_node_list) = [];

	        % partition state variables such that leaf nodes are assigned non-empty sets
	        % number of empty sets <= (number of non_leaf_nodes - 1)
	        % loop over number of empty sets
            
            num_empty = 0;
            nonempty_node_grp = node_list(1:(end-1));
            state_partitions = partitions(size(sys.X_DIMS_MASKED, 1), rr);
            assert(size(state_partitions{1}, 2)==size(nonempty_node_grp, 1), ...
		           'Number of state partitions must equal the number of non-empty nodes');
		    for spar = 1:1:size(state_partitions, 1)
                permutations_of_state_partitions = perms(state_partitions{spar, 1});
                for sper = 1:1:size(permutations_of_state_partitions, 1)
                    s = zeros(sys.U_DIMS, sys.X_DIMS);
                    for nne=1:1:size(nonempty_node_grp, 1)
                        state_indices = zeros(1, sys.X_DIMS);
                        state_indices(permutations_of_state_partitions{sper, nne}) = 1;
                        input_indices = rr_groupings{rrg}{nonempty_node_grp(nne)};
                        s(input_indices, :) = repmat(state_indices, [length(input_indices), 1]);
                    end

                    [c, ceq] = constraints(p, s);
                    assert(all(c <= 0) && all(abs(ceq) <= eps), 'Invalid decomposition!');
                    decomposition_count = decomposition_count + 1;
                    policy_decompositions{decomposition_count, 1} = p;
                    policy_decompositions{decomposition_count, 2} = s;
                    policy_decompositions{decomposition_count, 3} = computeLQRMeasure(sys, p, s);
                    policy_decompositions{decomposition_count, 4} = computeComplexityEstimates(sys, p, s);
                    policy_decompositions{decomposition_count, 5} = computeJointMeasure(sys, p, s);
                    fprintf('Decomposition : %d\n', decomposition_count);
                end
            end
            
	        for num_empty=1:1:(length(non_leaf_node_list) - 1)
	        	% pick 'empty' sets - nodes that will have no state assignment
	        	empty_node_groups = nchoosek(non_leaf_node_list(1:(end-1)), num_empty);
	        	if (~isempty(empty_node_groups))
		        	for ee = 1:1:size(empty_node_groups, 1)
		        		empty_node_grp = empty_node_groups(ee, :)';
		        		
		        		nonempty_node_grp = non_leaf_node_list(1:(end-1), :);
		        		nonempty_node_grp(any(nonempty_node_grp==empty_node_grp', 2)) = [];
                        if (isempty(nonempty_node_grp))
                            nonempty_node_grp = leaf_node_list;
                        else
                            nonempty_node_grp = cat(1, nonempty_node_grp, leaf_node_list);
                        end
                        
		        		state_partitions = partitions(size(sys.X_DIMS_MASKED, 1), rr - num_empty);
		        		assert(size(state_partitions{1}, 2)==size(nonempty_node_grp, 1), ...
		        			   'Number of state partitions must equal the number of non-empty nodes');
		        		if (isempty(state_partitions))
		        			continue;
		        		end
						for spar = 1:1:size(state_partitions, 1)
							permutations_of_state_partitions = perms(state_partitions{spar, 1});
							for sper = 1:1:size(permutations_of_state_partitions, 1)
								s = zeros(sys.U_DIMS, sys.X_DIMS);
								for nne=1:1:size(nonempty_node_grp, 1)
									state_indices = zeros(1, sys.X_DIMS);
                                    state_indices(permutations_of_state_partitions{sper, nne}) = 1;
                                    input_indices = rr_groupings{rrg}{nonempty_node_grp(nne)};
                                    s(input_indices, :) = repmat(state_indices, [length(input_indices), 1]);
								end

								[c, ceq] = constraints(p, s);
								assert(all(c <= 0) && all(abs(ceq) <= eps), 'Invalid decomposition!');
								decomposition_count = decomposition_count + 1;
								policy_decompositions{decomposition_count, 1} = p;
								policy_decompositions{decomposition_count, 2} = s;
								policy_decompositions{decomposition_count, 3} = computeLQRMeasure(sys, p, s);
								policy_decompositions{decomposition_count, 4} = computeComplexityEstimates(sys, p, s);
								policy_decompositions{decomposition_count, 5} = computeJointMeasure(sys, p, s);
                                fprintf('Decomposition : %d\n', decomposition_count);
							end
						end
		        	end
	        	end
	        end
		end
	end
end

