clear;
close all;
clc;

%% 

addpath('utils')
load('data/Biped2DSystemInputCoupled.mat');
NUM_PSEUDO_DIMS = size(sys.U_PSEUDO_DIMS, 1);

fun_casc = @(x) computeLQRMeasureCascaded(sys, reshape(x(1:NUM_PSEUDO_DIMS^2), NUM_PSEUDO_DIMS, NUM_PSEUDO_DIMS), ...
                                  reshape(x((NUM_PSEUDO_DIMS^2 + 1):(NUM_PSEUDO_DIMS^2 + NUM_PSEUDO_DIMS*sys.X_DIMS)), ...
                                          NUM_PSEUDO_DIMS, sys.X_DIMS));
nvars_casc = NUM_PSEUDO_DIMS^2 + sys.X_DIMS*NUM_PSEUDO_DIMS;
lb_casc = zeros(1, nvars_casc);
ub_casc = ones(1, nvars_casc);
A_casc = [];
b_casc = [];
Aeq_casc = [];
beq_casc = [];
nonlcon_casc = @(x) cascaded_constraints(reshape(x(1:NUM_PSEUDO_DIMS^2), NUM_PSEUDO_DIMS, NUM_PSEUDO_DIMS), ...
                                         reshape(x((NUM_PSEUDO_DIMS^2 + 1):(NUM_PSEUDO_DIMS^2 + NUM_PSEUDO_DIMS*sys.X_DIMS)), NUM_PSEUDO_DIMS, sys.X_DIMS));
IntCon_casc = linspace(1, nvars_casc, nvars_casc);
% options_casc = gaoptimset('OutputFcns', @debug_ga);
options_casc.Display = 'iter';
options_casc.PopulationSize = 500;
options_casc.CrossoverFraction = 0.9;
options_casc.EliteCount = 0.9*options_casc.PopulationSize;
% options_casc.PopulationType = 'bitstring';

[x_casc, err_lqr_casc, casc_exit_flag, casc_output, ...
 casc_population, casc_scores] = ga(fun_casc,nvars_casc,A_casc,b_casc,Aeq_casc,beq_casc,lb_casc,ub_casc,nonlcon_casc,IntCon_casc,options_casc);

% Extract decomposition
r_casc = round(reshape(x_casc(1:NUM_PSEUDO_DIMS^2), NUM_PSEUDO_DIMS, NUM_PSEUDO_DIMS));
s_casc = round(reshape(x_casc((NUM_PSEUDO_DIMS^2 + 1):(NUM_PSEUDO_DIMS^2 + NUM_PSEUDO_DIMS*sys.X_DIMS)), NUM_PSEUDO_DIMS, sys.X_DIMS));

%%

fun_dec = @(x) computeLQRMeasureDecoupled(sys, reshape(x, NUM_PSEUDO_DIMS, sys.X_DIMS));
nvars_dec = sys.X_DIMS*NUM_PSEUDO_DIMS;
lb_dec = zeros(1, nvars_dec);
ub_dec = ones(1, nvars_dec);
A_dec = [];
b_dec = [];
Aeq_dec = [];
beq_dec = [];
nonlcon_dec = @(x) decoupled_constraints(reshape(x, NUM_PSEUDO_DIMS, sys.X_DIMS));
IntCon_dec = linspace(1, nvars_dec, nvars_dec);

options_dec.Display = 'iter';
options_dec.PopulationSize = 500;
options_dec.CrossoverFraction = 0.9;
options_dec.EliteCount = 0.9*options_dec.PopulationSize;
% options_dec.PopulationType = 'bitstring';

[x_dec, err_lqr_dec, dec_exit_flag, dec_output, ...
 dec_population, dec_scores] = ga(fun_dec,nvars_dec,A_dec,b_dec,Aeq_dec,beq_dec,lb_dec,ub_dec,nonlcon_dec,IntCon_dec,options_dec);

% Extract decomposition
s_dec = reshape(x_dec(1:(NUM_PSEUDO_DIMS*sys.X_DIMS)), NUM_PSEUDO_DIMS, sys.X_DIMS);

%% Function

function [state, options, optchanged] = debug_ga(options, state, flag, interval)
    optchanged = false;
    [NUM_UNIQUE] = unique(state.Population, 'rows');
    load('data/Biped2DSystem.mat');
    num_invalid_samples = extract_decompositions_from_population(sys, state.Population);
end