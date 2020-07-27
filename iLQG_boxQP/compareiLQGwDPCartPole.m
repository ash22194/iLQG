clear;
% clc;
close all;

%%
addpath('cartpole');
iLQG_dir = 'data/';
iLQG_filenames = ["iLQGCartPoleDecomposed_diffR_difffurthercloserstarts_finestalpha_gamma_0.997mc=1,mp=1_1.mat"];
x_starts = [];

XJoint = [];
XCartFS = [];
XPoleTS = [];
XCartTS = [];
XPoleFS = [];
XCartFPoleTDec = [];
XCartTPoleFDec = [];
UJoint = [];
UCartFS = [];
UPoleTS = [];
UCartTS = [];
UPoleFS = [];
UCartFPoleTDec = [];
UCartTPoleFDec = [];

CJoint = [];
CPoleTS = [];
CCartFS = [];
CPoleFS = [];
CCartTS = [];

XCartFClose = [];
XPoleTClose = [];
XCartTClose = [];
XPoleFClose = [];
UCartFClose = [];
UPoleTClose = [];
UCartTClose = [];
UPoleFClose = [];
KCartFClose = [];
KPoleTClose = [];
KCartTClose = [];
KPoleFClose = [];

for ii=1:1:size(iLQG_filenames, 1)
    ilqg = load(strcat(iLQG_dir, iLQG_filenames(ii, 1)), 'XJoint', 'XCartFS', 'XPoleTS', 'XCartTS', 'XPoleFS', ...
                                      'XCartFPoleTDec', 'XCartTPoleFDec', ...
                                      'XCartFClose', 'XPoleTClose', 'XCartTClose', 'XPoleFClose', ...
                                      'UCartFClose', 'UPoleTClose', 'UCartTClose', 'UPoleFClose', ...
                                      'KCartFClose', 'KPoleTClose', 'KCartTClose', 'KPoleFClose', ...
                                      'UJoint', 'UCartFS', 'UPoleTS', 'UCartTS', 'UPoleFS', ...
                                      'UCartFPoleTDec', 'UCartTPoleFDec', ...
                                      'CJoint', 'CPoleTS', 'CCartFS', 'CPoleFS', 'CCartTS', ...
                                      'sys', 'x_starts');
    sys = ilqg.sys;
    x_starts(:,(end+1):(end+size(ilqg.x_starts, 2))) = ilqg.x_starts;
    
    XJoint = cat(3, XJoint, ilqg.XJoint);
    XCartFS = cat(3, XCartFS, ilqg.XCartFS);
    XPoleTS = cat(3, XPoleTS, ilqg.XPoleTS);
    XCartTS = cat(3, XCartTS, ilqg.XCartTS);
    XPoleFS = cat(3, XPoleFS, ilqg.XPoleFS);
    XCartFPoleTDec = cat(3, XCartFPoleTDec, ilqg.XCartFPoleTDec);
    XCartTPoleFDec = cat(3, XCartTPoleFDec, ilqg.XCartTPoleFDec);
    
    UJoint = cat(3, UJoint, ilqg.UJoint);
    UCartFS = cat(3, UCartFS, ilqg.UCartFS);
    UPoleTS = cat(3, UPoleTS, ilqg.UPoleTS);
    UCartTS = cat(3, UCartTS, ilqg.UCartTS);
    UPoleFS = cat(3, UPoleFS, ilqg.UPoleFS);
    UCartFPoleTDec = cat(3, UCartFPoleTDec, ilqg.UCartFPoleTDec);
    UCartTPoleFDec = cat(3, UCartTPoleFDec, ilqg.UCartTPoleFDec);
    
    CJoint  = cat(3, CJoint, ilqg.CJoint);
    CPoleTS = cat(3, CPoleTS, ilqg.CPoleTS);
    CCartFS = cat(3, CCartFS, ilqg.CCartFS);
    CPoleFS = cat(3, CPoleFS, ilqg.CPoleFS);
    CCartTS = cat(3, CCartTS, ilqg.CCartTS);
    
    XCartFClose = cat(3, XCartFClose, ilqg.XCartFClose);
    XPoleTClose = cat(3, XPoleTClose, ilqg.XPoleTClose);
    XCartTClose = cat(3, XCartTClose, ilqg.XCartTClose);
    XPoleFClose = cat(3, XPoleFClose, ilqg.XPoleFClose);
    
    UCartFClose = cat(3, UCartFClose, ilqg.UCartFClose);
    UPoleTClose = cat(3, UPoleTClose, ilqg.UPoleTClose);
    UCartTClose = cat(3, UCartTClose, ilqg.UCartTClose);
    UPoleFClose = cat(3, UPoleFClose, ilqg.UPoleFClose);
    
    KCartFClose = cat(4, KCartFClose, ilqg.KCartFClose);
    KPoleTClose = cat(4, KPoleTClose, ilqg.KPoleTClose);
    KCartTClose = cat(4, KCartTClose, ilqg.KCartTClose);
    KPoleFClose = cat(4, KPoleFClose, ilqg.KPoleFClose);
    
end

[x_starts, ia_x_starts, ~] = unique(x_starts', 'rows');
x_starts = x_starts';
% subset_starts = [1, 2, 5, 6, 7, 8, 11, 12, 19, 20, 23, 24, 25, 26, 29, 30];
% x_starts = x_starts(:, subset_starts);
% ia_x_starts = ia_x_starts(subset_starts);

XJoint = XJoint(:,:, ia_x_starts);
XCartFS = XCartFS(:,:, ia_x_starts);
XPoleTS = XPoleTS(:,:, ia_x_starts);
XCartTS = XCartTS(:,:, ia_x_starts);
XPoleFS = XPoleFS(:,:, ia_x_starts);
XCartFPoleTDec = XCartFPoleTDec(:,:, ia_x_starts);
XCartTPoleFDec = XCartTPoleFDec(:,:, ia_x_starts);

UJoint = UJoint(:,:, ia_x_starts);
UCartFS = UCartFS(:,:, ia_x_starts);
UPoleTS = UPoleTS(:,:, ia_x_starts);
UCartTS = UCartTS(:,:, ia_x_starts);
UPoleFS = UPoleFS(:,:, ia_x_starts);
UCartFPoleTDec = UCartFPoleTDec(:,:, ia_x_starts);
UCartTPoleFDec = UCartTPoleFDec(:,:, ia_x_starts);

CJoint  = CJoint(:,:, ia_x_starts);
CPoleTS = CPoleTS(:,:, ia_x_starts);
CCartFS = CCartFS(:,:, ia_x_starts);
CPoleFS = CPoleFS(:,:, ia_x_starts);
CCartTS = CCartTS(:,:, ia_x_starts);

XCartFClose = XCartFClose(:,:, ia_x_starts);
XPoleTClose = XPoleTClose(:,:, ia_x_starts);
XCartTClose = XCartTClose(:,:, ia_x_starts);
XPoleFClose = XPoleFClose(:,:, ia_x_starts);

UCartFClose = UCartFClose(:,:, ia_x_starts);
UPoleTClose = UPoleTClose(:,:, ia_x_starts);
UCartTClose = UCartTClose(:,:, ia_x_starts);
UPoleFClose = UPoleFClose(:,:, ia_x_starts);

KCartFClose = KCartFClose(:,:,:, ia_x_starts);
KPoleTClose = KPoleTClose(:,:,:, ia_x_starts);
KCartTClose = KCartTClose(:,:,:, ia_x_starts);
KPoleFClose = KPoleFClose(:,:,:, ia_x_starts);

DP_dir = 'data/';
DP_filename = 'mc=1,mp=1_densegrid_gamma_0.997.mat';

load(strcat(DP_dir, DP_filename), 'V_joint', 'policy_joint', ...
                                  'V_t_pole_second', 'f_cart_first_', 't_pole_second', ...
                                  'V_f_cart_second', 't_pole_first_', 'f_cart_second', ...
                                  'V_f_pole_second', 't_cart_first_', 'f_pole_second', ...
                                  'V_t_cart_second', 'f_pole_first_', 't_cart_second', ...
                                  'V_f_cart_t_pole', 'V_t_cart_f_pole', ...
                                  'grid_x', 'grid_x_dot', 'grid_xP', 'grid_xP_dot', ...
                                  'x_limits', 'x_dot_limits', 'xP_limits', 'xP_dot_limits');

%% Estimate DDP values
J_DDP = nan(size(x_starts,2), 7);
J_DDP_final = nan(size(x_starts,2), 5);

DDP_finalerr = nan(size(x_starts,2), 7);
convergence_tol = 5e-2;

disp('DDP');
for jj=1:1:size(x_starts, 2)
    
    J_DDP(jj, 1)         = J_CartPole(sys, XJoint(:,:, jj), UJoint(:,:, jj));
    J_DDP_final(jj, 1)   = sum(CJoint(1,:, jj), 2);
    lastState            = XJoint(:, end, jj);
    lastState(3)         = mod(lastState(3), 2*pi);
    DDP_finalerr(jj, 1)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [2];
    sys.U_DIMS_FIXED     = [1];
    cost                 = l_CartPoleSecond(sys, XPoleTS(:,:, jj), UPoleTS(:,:, jj), ...
                                              UCartFClose(:,:, jj), KCartFClose(:,:,:, jj), XCartFClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 2)         = sum(cost.*discount*sys.dt, 2);
    J_DDP_final(jj, 2)   = sum(CPoleTS(1,:, jj), 2);
    lastState            = XPoleTS(:, end, jj);
    lastState(3)         = mod(lastState(3), 2*pi);
    DDP_finalerr(jj, 2)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [1];
    sys.U_DIMS_FIXED     = [2];
    cost                 = l_CartPoleSecond(sys, XCartFS(:,:, jj), UCartFS(:,:, jj), ...
                                              UPoleTClose(:,:, jj), KPoleTClose(:,:,:, jj), XPoleTClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 3)         = sum(cost.*discount*sys.dt, 2);
    J_DDP_final(jj, 3)   = sum(CCartFS(1,:, jj), 2);
    lastState            = XCartFS(:, end, jj);
    lastState(3)         = mod(lastState(3), 2*pi);
    DDP_finalerr(jj, 3)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [1];
    sys.U_DIMS_FIXED     = [2];
    cost                 = l_CartPoleSecond(sys, XPoleFS(:,:, jj), UPoleFS(:,:, jj), ...
                                              UCartTClose(:,:, jj), KCartTClose(:,:,:, jj), XCartTClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 4)         = sum(cost.*discount*sys.dt, 2);
    J_DDP_final(jj, 4)   = sum(CPoleFS(1,:, jj), 2);
    lastState            = XPoleFS(:, end, jj);
    lastState(3)         = mod(lastState(3), 2*pi);
    DDP_finalerr(jj, 4)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [2];
    sys.U_DIMS_FIXED     = [1];
    cost                 = l_CartPoleSecond(sys, XCartTS(:,:, jj), UCartTS(:,:, jj), ...
                                              UPoleFClose(:,:, jj), KPoleFClose(:,:,:, jj), XPoleFClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 5)         = sum(cost.*discount*sys.dt, 2);
    J_DDP_final(jj, 5)   = sum(CCartTS(1,:, jj), 2);
    lastState            = XCartTS(:, end, jj);
    lastState(3)         = mod(lastState(3), 2*pi);
    DDP_finalerr(jj, 5)  = (norm(lastState - sys.goal));
    
    J_DDP(jj, 6)         = J_CartPole(sys, XCartFPoleTDec(:,:, jj), UCartFPoleTDec(:,:, jj));
    lastState            = XCartFPoleTDec(:, end, jj);
    lastState(3)         = mod(lastState(3), 2*pi);
    DDP_finalerr(jj, 6)  = (norm(lastState - sys.goal));
    
    J_DDP(jj, 7)         = J_CartPole(sys, XCartTPoleFDec(:,:, jj), UCartTPoleFDec(:,:, jj));
    lastState            = XCartTPoleFDec(:, end, jj);
    lastState(3)         = mod(lastState(3), 2*pi);
    DDP_finalerr(jj, 7)  = (norm(lastState - sys.goal));
    jj
end
DDP_converged = (DDP_finalerr < convergence_tol);

%% DP
disp('DP');
sys.U_DIMS_FREE = [1;2];
sys.U_DIMS_FIXED = [];
sys.X_DIMS_FREE = [1;2;3;4];
sys.X_DIMS_FIXED = [];
limits = [x_limits; x_dot_limits; xP_limits; xP_dot_limits];

J_DP_rollout = nan(size(x_starts, 2), 7);
J_DP_longhorz_rollout = nan(size(x_starts, 2), 7);
DP_finalerr = nan(size(x_starts, 2), 7);
DP_longhorz_finalerr = nan(size(x_starts, 2), 7);
T = 5;
NUM_STEPS = round(T / sys.dt);
T_longhorz = 10;
NUM_STEPS_LONGHORZ = round(T_longhorz / sys.dt);

disp('Joint');
[trajXJoint, trajUJoint, J_DP_rollout(:,1)] = rolloutDPTrajOde(policy_joint, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
trajXJointErr = trajXJoint(:, end, :) - sys.goal;
DP_finalerr(:, 1) = reshape(vecnorm(trajXJointErr), size(x_starts,2), 1);

[trajXJoint, trajUJoint, J_DP_longhorz_rollout(:,1)] = rolloutDPTrajOde(policy_joint, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
trajXJointErr = trajXJoint(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 1) = reshape(vecnorm(trajXJointErr), size(x_starts,2), 1);

disp('Dec 1');
policyCartFPoleT(:,:,:,:,1) = f_cart_first_;
policyCartFPoleT(:,:,:,:,2) = t_pole_second;

[trajXCFPT, trajUCFPT, J_DP_rollout(:,2)] = rolloutDPTrajOde(policyCartFPoleT, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
trajXCFPTErr = trajXCFPT(:, end, :) - sys.goal;
DP_finalerr(:, 2) = reshape(vecnorm(trajXCFPTErr), size(x_starts,2), 1);

[trajXCFPT, trajUCFPT, J_DP_longhorz_rollout(:,2)] = rolloutDPTrajOde(policyCartFPoleT, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
trajXCFPTErr = trajXCFPT(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 2) = reshape(vecnorm(trajXCFPTErr), size(x_starts,2), 1);

disp('Dec 2');
policyPoleTCartF(:,:,:,:,2) = t_pole_first_;
policyPoleTCartF(:,:,:,:,1) = f_cart_second;

[trajXPTCF, trajUPTCF, J_DP_rollout(:,3)] = rolloutDPTrajOde(policyPoleTCartF, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
trajXPTCFErr = trajXPTCF(:, end, :) - sys.goal;
DP_finalerr(:, 3) = reshape(vecnorm(trajXPTCFErr), size(x_starts,2), 1);

[trajXPTCF, trajUPTCF, J_DP_longhorz_rollout(:,3)] = rolloutDPTrajOde(policyPoleTCartF, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
trajXPTCFErr = trajXPTCF(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 3) = reshape(vecnorm(trajXPTCFErr), size(x_starts,2), 1);

disp('Dec 3');
policyCartTPoleF(:,:,:,:,2) = t_cart_first_;
policyCartTPoleF(:,:,:,:,1) = f_pole_second;

[trajXCTPF, trajUCTPF, J_DP_rollout(:,4)] = rolloutDPTrajOde(policyCartTPoleF, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
trajXCTPFErr = trajXCTPF(:, end, :) - sys.goal;
DP_finalerr(:, 4) = reshape(vecnorm(trajXCTPFErr), size(x_starts,2), 1);

[trajXCTPF, trajUCTPF, J_DP_longhorz_rollout(:,4)] = rolloutDPTrajOde(policyCartTPoleF, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
trajXCTPFErr = trajXCTPF(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 4) = reshape(vecnorm(trajXCTPFErr), size(x_starts,2), 1);

disp('Dec 4');
policyPoleFCartT(:,:,:,:,1) = f_pole_first_;
policyPoleFCartT(:,:,:,:,2) = t_cart_second;
[trajXPFCT, trajUPFCT, J_DP_rollout(:,5)] = rolloutDPTrajOde(policyPoleFCartT, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
trajXPFCTErr = trajXPFCT(:, end, :) - sys.goal;
DP_finalerr(:, 5) = reshape(vecnorm(trajXPFCTErr), size(x_starts,2), 1);

[trajXPFCT, trajUPFCT, J_DP_longhorz_rollout(:,5)] = rolloutDPTrajOde(policyPoleFCartT, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
trajXPFCTErr = trajXPFCT(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 5) = reshape(vecnorm(trajXPFCTErr), size(x_starts,2), 1);

disp('Dec 5');
policyCartFPoleTDec(:,:,:,:,1) = f_cart_first_;
policyCartFPoleTDec(:,:,:,:,2) = t_pole_first_;

[trajXCFPTD, trajUCFPTD, J_DP_rollout(:,6)] = rolloutDPTrajOde(policyCartFPoleTDec, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
trajXCFPTDErr = trajXCFPTD(:, end, :) - sys.goal;
DP_finalerr(:, 6) = reshape(vecnorm(trajXCFPTDErr), size(x_starts,2), 1);

[trajXCFPTD, trajUCFPTD, J_DP_longhorz_rollout(:,6)] = rolloutDPTrajOde(policyCartFPoleTDec, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
trajXCFPTDErr = trajXCFPTD(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 6) = reshape(vecnorm(trajXCFPTDErr), size(x_starts,2), 1);

disp('Dec 6');
policyCartTPoleFDec(:,:,:,:,1) = f_pole_first_;
policyCartTPoleFDec(:,:,:,:,2) = t_cart_first_;

[trajXCTPFD, trajUCTPFD, J_DP_rollout(:,7)] = rolloutDPTrajOde(policyCartTPoleFDec, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
trajXCTPFDErr = trajXCTPFD(:, end, :) - sys.goal;
DP_finalerr(:, 7) = reshape(vecnorm(trajXCTPFDErr), size(x_starts,2), 1);

[trajXCTPFD, trajUCTPFD, J_DP_longhorz_rollout(:,7)] = rolloutDPTrajOde(policyCartTPoleFDec, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
trajXCTPFDErr = trajXCTPFD(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 7) = reshape(vecnorm(trajXCTPFDErr), size(x_starts,2), 1);

DP_converged = DP_finalerr < convergence_tol;
DP_longhorz_converged = DP_longhorz_finalerr < convergence_tol;

%% DP values for starts
disp('DP');
J_DP = nan(size(x_starts, 2), 7);
J_DP(:, 1) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, V_joint, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)');
                 
J_DP(:, 2) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, V_t_pole_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)');

J_DP(:, 3) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, V_f_cart_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)');

J_DP(:, 4) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, V_f_pole_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)');

J_DP(:, 5) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, V_t_cart_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)');

J_DP(:, 6) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, V_f_cart_t_pole, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)');

J_DP(:, 7) = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, V_t_cart_f_pole, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)');

%% Ordering

% DDP
DDP_err = sum(abs(J_DDP(:,1) - J_DDP(:,2:end)));
[~,DDP_order] = sort(DDP_err);

% DP
DP_err = sum(abs(J_DP(:,1) - J_DP(:,2:end)));
[~,DP_order] = sort(DP_err);

% DP rollout
DP_rollout_err = sum(abs(J_DP_rollout(:,1) - J_DP_rollout(:,2:end)));
[~,DP_rollout_order] = sort(DP_rollout_err);
DP_longhorz_rollout_err = sum(abs(J_DP_longhorz_rollout(:,1) - J_DP_longhorz_rollout(:,2:end)));
[~,DP_longhorz_rollout_order] = sort(DP_longhorz_rollout_err);

%% DP Value function comparison within state bounds

state_bounds = [-0.5, 0.5;
                -1, 1;
                2*pi/3, 4*pi/3;
                -1, 1];

valid_range = ((grid_x >= state_bounds(1,1)) & (grid_x <= state_bounds(1,2)) ...
                & (grid_x_dot >= state_bounds(2,1)) & (grid_x_dot <= state_bounds(2,2)) ...
                & (grid_xP >= state_bounds(3,1)) & (grid_xP <= state_bounds(3,2)) ...
                & (grid_xP_dot >= state_bounds(4,1)) & (grid_xP_dot <= state_bounds(4,2)));

DP_avg_err = [sum(abs(V_joint(valid_range) - V_t_pole_second(valid_range)), 'all'), ... 
              sum(abs(V_joint(valid_range) - V_f_cart_second(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_f_pole_second(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_t_cart_second(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_f_cart_t_pole(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_t_cart_f_pole(valid_range)), 'all')];

[~, DP_avg_order] = sort(DP_avg_err);

%% Functions
function [trajX, trajU, J, J_WrapAround] = rolloutDPTraj(policy, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T)
    
    NUM_CTRL = round(T / sys.dt);
    trajX = nan(4, NUM_CTRL+1, size(x_starts, 2));
    trajU = nan(2, NUM_CTRL, size(x_starts, 2));
    J = zeros(size(x_starts, 2), 1);
    J_WrapAround = zeros(size(x_starts, 2), 1);
    P1 = griddedInterpolant(grid_x, grid_x_dot, grid_xP, grid_xP_dot, policy(:,:,:,:,1));
    P2 = griddedInterpolant(grid_x, grid_x_dot, grid_xP, grid_xP_dot, policy(:,:,:,:,2));
    
    for jj = 1:1:size(x_starts, 2)
        trajX(:, 1, jj) = x_starts(: ,jj);
        for ii=1:1:NUM_CTRL
            x_ = [min(limits(1, 2), max(limits(1, 1), trajX(1, ii, jj)));
                  min(limits(2, 2), max(limits(2, 1), trajX(2, ii, jj)));
                  min(limits(3, 2), max(limits(3, 1), trajX(3, :, jj)));
%                   mod(trajX(3, :, jj), 2*pi);
                  min(limits(4, 2), max(limits(4, 1), trajX(4, ii, jj)))];
            
            trajU(1, ii, jj) = P1(x_(1),  x_(2),  x_(3),  x_(4));
            trajU(2, ii, jj) = P2(x_(1),  x_(2),  x_(3),  x_(4));
            [trajX(:, ii+1, jj),~] = sys.dynamics_discrete(trajX(:, ii, jj), trajU(:, ii, jj));
        end
        J_WrapAround(jj, 1) = J_CartPole_WrapAround(sys, trajX(:,:, jj), trajU(:,:, jj));
        J(jj, 1) = J_CartPole(sys, trajX(:,:, jj), trajU(:,:, jj));
    end
end

function [trajX, trajU, J] = rolloutDPTrajOde(policy, x_starts, sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T)
    
    NUM_CTRL = round(T / sys.dt);
    tspan = linspace(0, T, NUM_CTRL+1);
    opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
    trajX = nan(4, NUM_CTRL+1, size(x_starts, 2));
    trajU = nan(2, NUM_CTRL, size(x_starts, 2));
    J = zeros(size(x_starts, 2), 1);
    P1 = griddedInterpolant(grid_x, grid_x_dot, grid_xP, grid_xP_dot, policy(:,:,:,:,1));
    P2 = griddedInterpolant(grid_x, grid_x_dot, grid_xP, grid_xP_dot, policy(:,:,:,:,2));
    for jj = 1:1:size(x_starts, 2)
        [~, X] = ode113(@(t,y) cartpole_dyn_gridbased(t, y, P1, P2, limits, sys), tspan, x_starts(:, jj), opts);
        trajX(:,:, jj) = X';
        X_ = [min(limits(1, 2), max(limits(1, 1), trajX(1, :, jj)));
              min(limits(2, 2), max(limits(2, 1), trajX(2, :, jj)));
              min(limits(3, 2), max(limits(3, 1), trajX(3, :, jj)));
              min(limits(4, 2), max(limits(4, 1), trajX(4, :, jj)))];
        
        trajU(1, :, jj) = P1(X_(1,1:(end-1)), X_(2,1:(end-1)), X_(3,1:(end-1)), X_(4,1:(end-1)));
        trajU(2, :, jj) = P2(X_(1,1:(end-1)), X_(2,1:(end-1)), X_(3,1:(end-1)), X_(4,1:(end-1)));
        J(jj, 1) = J_CartPole(sys, trajX(:,:, jj), trajU(:,:, jj));
        disp(strcat('start : ', num2str(jj)));
    end
    
end