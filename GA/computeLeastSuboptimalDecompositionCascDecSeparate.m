clear;
close all;
clc;

%% 

addpath('utils');
load('data/Biped2DSystem.mat');
NUM_PSEUDO_DIMS = size(sys.U_PSEUDO_DIMS, 1);

fun_casc = @(x) computeLQRMeasureCascaded(sys, reshape(x(1:NUM_PSEUDO_DIMS^2), NUM_PSEUDO_DIMS, NUM_PSEUDO_DIMS), ...
                                  reshape(x((NUM_PSEUDO_DIMS^2 + 1):(NUM_PSEUDO_DIMS^2 + NUM_PSEUDO_DIMS*sys.X_DIMS)), ...
                                          NUM_PSEUDO_DIMS, sys.X_DIMS));
nvars_casc = NUM_PSEUDO_DIMS^2 + sys.X_DIMS*NUM_PSEUDO_DIMS;
lb_casc = zeros(1, nvars_casc);
ub_casc = ones(1, nvars_casc);
Aeq_casc = [];
beq_casc = [];
nonlcon_casc = [];
IntCon_casc = linspace(1,nvars_casc,nvars_casc);
% options_casc = gaoptimset('OutputFcns', @debug_ga);
options_casc.Display = 'iter';
options_casc.PopulationSize = 500;
options_casc.CrossoverFraction = 0.9;
options_casc.EliteCount = 0.9*options_casc.PopulationSize;
% options_casc.PopulationType = 'bitstring';

nconst_casc = 4*(NUM_PSEUDO_DIMS) + NUM_PSEUDO_DIMS + 2*sys.X_DIMS;
A_casc = zeros(nconst_casc, nvars_casc);
b_casc = zeros(nconst_casc, 1);

% Unique order constraint
for ii=1:1:NUM_PSEUDO_DIMS
    % One rank is assigned to only one input
    A_casc(ii, ((ii-1)*NUM_PSEUDO_DIMS + 1):(ii*NUM_PSEUDO_DIMS)) = ones(1, NUM_PSEUDO_DIMS);
    b_casc(ii) = 1;
    A_casc(NUM_PSEUDO_DIMS + ii, ((ii-1)*NUM_PSEUDO_DIMS + 1):(ii*NUM_PSEUDO_DIMS)) = -ones(1, NUM_PSEUDO_DIMS);
    b_casc(NUM_PSEUDO_DIMS + ii) = -1;
    
    % Each input is assigned a rank
    A_casc(2*NUM_PSEUDO_DIMS + ii, ii + NUM_PSEUDO_DIMS*linspace(0,NUM_PSEUDO_DIMS-1,NUM_PSEUDO_DIMS)) = ones(1, NUM_PSEUDO_DIMS);
    b_casc(2*NUM_PSEUDO_DIMS + ii) = 1;
    A_casc(3*NUM_PSEUDO_DIMS + ii, ii + NUM_PSEUDO_DIMS*linspace(0,NUM_PSEUDO_DIMS-1,NUM_PSEUDO_DIMS)) = -ones(1, NUM_PSEUDO_DIMS);
    b_casc(3*NUM_PSEUDO_DIMS + ii) = -1;
    
    % Each input must have at least one state variable
    A_casc(4*NUM_PSEUDO_DIMS + ii, NUM_PSEUDO_DIMS^2 + ii + NUM_PSEUDO_DIMS*linspace(0,sys.X_DIMS-1,sys.X_DIMS)) = -ones(1, sys.X_DIMS);
end

for jj=1:1:sys.X_DIMS
    % Each state variable must be assigned to ony one input
    A_casc(5*NUM_PSEUDO_DIMS + jj, ...
      (NUM_PSEUDO_DIMS^2 + ((jj-1)*NUM_PSEUDO_DIMS + 1)):(NUM_PSEUDO_DIMS^2 + jj*NUM_PSEUDO_DIMS)) = ones(1, NUM_PSEUDO_DIMS);
    b_casc(5*NUM_PSEUDO_DIMS + jj) = 1;
    A_casc(5*NUM_PSEUDO_DIMS + sys.X_DIMS + jj, ...
      (NUM_PSEUDO_DIMS^2 + ((jj-1)*NUM_PSEUDO_DIMS + 1)):(NUM_PSEUDO_DIMS^2 + jj*NUM_PSEUDO_DIMS)) = -ones(1, NUM_PSEUDO_DIMS);
    b_casc(5*NUM_PSEUDO_DIMS + sys.X_DIMS + jj) = -1;
end

[x_casc, err_lqr_casc, casc_exit_flag, casc_output, ...
 casc_population, casc_scores] = ga(fun_casc,nvars_casc,A_casc,b_casc,Aeq_casc,beq_casc,lb_casc,ub_casc,nonlcon_casc,IntCon_casc,options_casc);

% Extract decomposition
r_casc = reshape(x_casc(1:NUM_PSEUDO_DIMS^2), NUM_PSEUDO_DIMS, NUM_PSEUDO_DIMS);
s_casc = reshape(x_casc((NUM_PSEUDO_DIMS^2 + 1):(NUM_PSEUDO_DIMS^2 + NUM_PSEUDO_DIMS*sys.X_DIMS)), NUM_PSEUDO_DIMS, sys.X_DIMS);

%%

fun_dec = @(x) computeLQRMeasureDecoupled(sys, reshape(x, NUM_PSEUDO_DIMS, sys.X_DIMS));
nvars_dec = sys.X_DIMS*NUM_PSEUDO_DIMS;
lb_dec = zeros(1, nvars_dec);
ub_dec = ones(1, nvars_dec);
Aeq_dec = [];
beq_dec = [];
nonlcon_dec = [];
IntCon_dec = linspace(1,nvars_dec,nvars_dec);
options_dec.Display = 'iter';
options_dec.PopulationSize = 500;
options_dec.CrossoverFraction = 0.9;
options_dec.EliteCount = 0.9*options_dec.PopulationSize;
% options_dec.PopulationType = 'bitstring';
nconst_dec = NUM_PSEUDO_DIMS + 2*sys.X_DIMS;
A_dec = zeros(nconst_dec, nvars_dec);
b_dec = zeros(nconst_dec, 1);

% Each input must have at least one state variable
for ii=1:1:NUM_PSEUDO_DIMS
    A_dec(ii, ii + NUM_PSEUDO_DIMS*linspace(0,sys.X_DIMS-1,sys.X_DIMS)) = -ones(1, sys.X_DIMS);
    b_dec(ii) = -1;
end

for jj=1:1:sys.X_DIMS
    % Each state variable must be assigned to ony one input
    A_dec(NUM_PSEUDO_DIMS + jj, ((jj-1)*NUM_PSEUDO_DIMS + 1):(jj*NUM_PSEUDO_DIMS)) = ones(1, NUM_PSEUDO_DIMS);
    b_dec(NUM_PSEUDO_DIMS + jj) = 1;
    A_dec(NUM_PSEUDO_DIMS + sys.X_DIMS + jj, ...
         ((jj-1)*NUM_PSEUDO_DIMS + 1):(jj*NUM_PSEUDO_DIMS)) = -ones(1, NUM_PSEUDO_DIMS);
    b_dec(NUM_PSEUDO_DIMS + sys.X_DIMS + jj) = -1;
end

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