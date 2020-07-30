clear;
close all;
clc;

%% 

load('data/Biped2DSystem.mat');
U_DIMS = sys.U_DIMS;

% fun_casc = @(x) computeLQRMeasureCascadedwPseudoInputs(sys, reshape(x(1:U_DIMS^2), U_DIMS, U_DIMS), ...
%                                   reshape(x((U_DIMS^2 + 1):(U_DIMS^2 + U_DIMS*sys.X_DIMS)), U_DIMS, sys.X_DIMS));
% nvars_casc = U_DIMS^2 + sys.X_DIMS*U_DIMS;
fun_casc = @(x) computeLQRMeasureCascadedwRealRanks(sys, reshape(x(1:U_DIMS), U_DIMS, 1), ...
                                  reshape(x((U_DIMS + 1):(U_DIMS + U_DIMS*sys.X_DIMS)), U_DIMS, sys.X_DIMS));
nvars_casc = U_DIMS + sys.X_DIMS*U_DIMS;
lb_casc = zeros(1, nvars_casc);
lb_casc(1:U_DIMS) = 1;
ub_casc = ones(1, nvars_casc);
ub_casc(1:U_DIMS) = U_DIMS;
A_casc = [];
b_casc = [];
Aeq_casc = [];
beq_casc = [];
% nonlcon_casc = @(x) cascaded_constraints_wpseudo_inputs(reshape(x(1:U_DIMS^2), U_DIMS, U_DIMS), ...
%                                   reshape(x((U_DIMS^2 + 1):(U_DIMS^2 + U_DIMS*sys.X_DIMS)), U_DIMS, sys.X_DIMS));
nonlcon_casc = @(x) cascaded_constraints_wreal_ranks(reshape(x(1:U_DIMS), U_DIMS, 1), ...
                                  reshape(x((U_DIMS + 1):(U_DIMS + U_DIMS*sys.X_DIMS)), U_DIMS, sys.X_DIMS));
IntCon_casc = []; % linspace(1,nvars_casc,nvars_casc);

% options_casc = gaoptimset('OutputFcns', @debug_ga);
options_casc.Display = 'iter';
options_casc.PopulationSize = 500;
options_casc.CrossoverFraction = 0.9;
options_casc.EliteCount = 0.9*options_casc.PopulationSize;
options_casc.CreationFcn = @(nvars, fitness_fcn, options) generate_population(sys, options.PopulationSize);
% options_casc.PopulationType = 'bitstring';

[x_casc, err_lqr_casc, casc_exit_flag, casc_output, ...
 casc_population, casc_scores] = ga(fun_casc,nvars_casc,A_casc,b_casc,Aeq_casc,beq_casc,lb_casc,ub_casc,nonlcon_casc,IntCon_casc,options_casc);
% x0 = generate_population(sys, 1);
% [x_casc, err_lqr_casc, casc_exit_flag, casc_output] = patternsearch(fun_casc,x0,A_casc,b_casc,Aeq_casc,beq_casc,lb_casc,ub_casc,nonlcon_casc,options_casc);

% Extract decomposition
% r_casc = reshape(x_casc(1:U_DIMS^2), U_DIMS, U_DIMS);
% s_casc = reshape(x_casc((U_DIMS^2 + 1):(U_DIMS^2 + U_DIMS*sys.X_DIMS)), U_DIMS, sys.X_DIMS);
r_casc = round(reshape(x_casc(1:U_DIMS), U_DIMS, 1));
s_casc = round(reshape(x_casc((U_DIMS + 1):(U_DIMS + U_DIMS*sys.X_DIMS)), U_DIMS, sys.X_DIMS));

%%

fun_dec = @(x) computeLQRMeasureDecoupled(sys, reshape(x, U_DIMS, sys.X_DIMS));
nvars_dec = sys.X_DIMS*U_DIMS;
lb_dec = zeros(1, nvars_dec);
ub_dec = ones(1, nvars_dec);
A_dec = [];
b_dec = [];
Aeq_dec = [];
beq_dec = [];
nonlcon_dec = @(x) decoupled_constraints_wpseudo_inputs(reshape(x, U_DIMS, sys.X_DIMS));
IntCon_dec = linspace(1,nvars_dec,nvars_dec);

options_dec.Display = 'iter';
options_dec.PopulationSize = 1000;
options_dec.CrossoverFraction = 0.9;
options_dec.EliteCount = 0.9*options_dec.PopulationSize;
% options_dec.PopulationType = 'bitstring';

[x_dec, err_lqr_dec, dec_exit_flag, dec_output, ...
 dec_population, dec_scores] = ga(fun_dec,nvars_dec,A_dec,b_dec,Aeq_dec,beq_dec,lb_dec,ub_dec,nonlcon_dec,IntCon_dec,options_dec);

% Extract decomposition
s_dec = reshape(x_dec(1:(U_DIMS*sys.X_DIMS)), U_DIMS, sys.X_DIMS);

%% Function

function [state, options, optchanged] = debug_ga(options, state, flag, interval)
    optchanged = false;
    [NUM_UNIQUE] = unique(state.Population, 'rows');
    load('data/Biped2DSystem.mat');
    num_invalid_samples = extract_decompositions_from_population(sys, state.Population);
end

function population = generate_population(sys, n)
    
    r = randi([1,sys.U_DIMS], sys.U_DIMS, n);
    s = zeros(sys.U_DIMS, sys.X_DIMS, n);
    for ii=1:1:n
        rC = unique(r(:,ii));
        state = linspace(1, sys.X_DIMS, sys.X_DIMS);
        for jj=1:1:length(rC)
            k = randi([1, length(state)-(length(rC) - jj)], 1);
            sub_state = nchoosek(linspace(1,length(state),length(state)), k);
            sub_state = sub_state(randi([1, size(sub_state,1)],1),:);
            sub_state_ = state(sub_state);
            state(sub_state) = [];

            s(r(:,ii)==rC(jj), sub_state_, ii) = 1;
        end
        s(r(:,ii)==rC(jj), state, ii) = 1;
    end
    
    population = [r', reshape(permute(s,[3,1,2]), n, sys.U_DIMS*sys.X_DIMS)];
end