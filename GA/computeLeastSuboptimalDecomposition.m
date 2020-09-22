clear;
close all;
clc;

%% 

system_name = 'manipulator2dof';
load(strcat('data/', system_name, 'System.mat'));
sys.xu = [sys.x; sys.u];

fun = @(x) computeLQRMeasure(sys, reshape(x(1:2*sys.U_DIMS), sys.U_DIMS, 2), ...
                                  reshape(x((2*sys.U_DIMS + 1):(2*sys.U_DIMS + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS));
nvars = 2*sys.U_DIMS + sys.X_DIMS*sys.U_DIMS;
lb = [zeros(1, sys.U_DIMS), ones(1, sys.U_DIMS), zeros(1, sys.U_DIMS*sys.X_DIMS)];
ub = [sys.U_DIMS*ones(1, sys.U_DIMS*2), ones(1, sys.U_DIMS*sys.X_DIMS)];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = @(x) constraints(reshape(x(1:2*sys.U_DIMS), sys.U_DIMS, 2), ...
                           reshape(x((2*sys.U_DIMS + 1):(2*sys.U_DIMS + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS));
IntCon = [];
options.Display = 'iter';
options.PopulationSize = 200;
options.CrossoverFraction = 0.6;
options.EliteCount = 0.2*options.PopulationSize;
options.InitialPopulation = generate_population(sys, options.PopulationSize);
for ii=1:1:options.PopulationSize
    ii
    p = reshape(options.InitialPopulation(ii, 1:2*sys.U_DIMS), sys.U_DIMS,2);
    s = reshape(options.InitialPopulation(ii, (1+2*sys.U_DIMS):end), sys.U_DIMS,sys.X_DIMS);
    [c, c_eq] = constraints(p, s);
    if (any(c_eq~=0) || any(c > 0))
        disp('Invalid sample');
    end
end
options.CreationFcn = @(nvars, fitness_fcn, options) ...
                        generate_population(sys, options.PopulationSize);
options.CrossoverFcn = @(parents, options, nvars, fitness_fcn, unused, population) ...
                         crossoverfunction(sys, parents, options, nvars, fitness_fcn, unused, population);
options.MutationFcn = @(parents, options, nvars, fitness_fcn, state, score, population) ...
                         mutationfunction(sys, parents, options, nvars, fitness_fcn, state, score, population);
profile on;
[x, err_lqr, exitflag, output, population, scores] = ga(fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,IntCon,options);
profile off;
profile viewer;

% Extract decomposition
p = reshape(x(1:2*sys.U_DIMS), sys.U_DIMS, 2);
s = reshape(x((2*sys.U_DIMS + 1):(2*sys.U_DIMS + sys.U_DIMS*sys.X_DIMS)), sys.U_DIMS, sys.X_DIMS);

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