clear;
% clc;
close all;

%%
addpath('biped2d');
iLQG_dir = 'data/';
iLQG_filenames = ["iLQGBiped2DDecomposed_newerthetastarts_corrlin_all8.mat";
                  ];
x_starts = [];

XJoint = [];
XCOMFS = [];
XTorsoTS = [];
XCOMTS = [];
XTorsoFS = [];
XTSFull = [];
XFSFull = [];
XCOMFTorsoTDec = [];
XCOMTTorsoFDec = [];

UJoint = [];
UCOMFS = [];
UTorsoTS = [];
UCOMTS = [];
UTorsoFS = [];
UTSFull = [];
UFSFull = [];
UCOMFTorsoTDec = [];
UCOMTTorsoFDec = [];

% CJoint = [];
% CTorsoTS = [];
% CCOMFS = [];
% CTorsoFS = [];
% CCOMTS = [];

XCOMFClose = [];
XTorsoTClose = [];
XCOMTClose = [];
XTorsoFClose = [];
XFFullClose = [];
XTFullClose = [];

UCOMFClose = [];
UTorsoTClose = [];
UCOMTClose = [];
UTorsoFClose = [];
UFFullClose = [];
UTFullClose = [];

KCOMFClose = [];
KTorsoTClose = [];
KCOMTClose = [];
KTorsoFClose = [];
KFFullClose = [];
KTFullClose = [];

for ii=1:1:size(iLQG_filenames, 1)
    ilqg = load(strcat(iLQG_dir, iLQG_filenames(ii, 1)), 'XJoint', 'XCOMFS', 'XTorsoTS', 'XCOMTS', 'XTorsoFS', 'XTSFull', 'XFSFull', ...
                                      'XCOMFTorsoTDec', 'XCOMTTorsoFDec', ...
                                      'XCOMFClose', 'XTorsoTClose', 'XCOMTClose', 'XTorsoFClose', 'XFFullClose', 'XTFullClose', ...
                                      'UCOMFClose', 'UTorsoTClose', 'UCOMTClose', 'UTorsoFClose', 'UFFullClose', 'UTFullClose', ...
                                      'KCOMFClose', 'KTorsoTClose', 'KCOMTClose', 'KTorsoFClose', 'KFFullClose', 'KTFullClose', ...
                                      'UJoint', 'UCOMFS', 'UTorsoTS', 'UCOMTS', 'UTorsoFS', 'UTSFull', 'UFSFull', ...
                                      'UCOMFTorsoTDec', 'UCOMTTorsoFDec', ...
                                      'CJoint', 'CTorsoTS', 'CCOMFS', 'CTorsoFS', 'CCOMTS', 'CTSFull', 'CFSFull', ...
                                      'sys', 'x_starts');
    sys = ilqg.sys;
    x_starts(:,(end+1):(end+size(ilqg.x_starts, 2))) = ilqg.x_starts;
    
    XJoint = cat(3, XJoint, ilqg.XJoint);
    XCOMFS = cat(3, XCOMFS, ilqg.XCOMFS);
    XTorsoTS = cat(3, XTorsoTS, ilqg.XTorsoTS);
    XCOMTS = cat(3, XCOMTS, ilqg.XCOMTS);
    XTorsoFS = cat(3, XTorsoFS, ilqg.XTorsoFS);
    XTSFull = cat(3, XTSFull, ilqg.XTSFull);
    XFSFull = cat(3, XFSFull, ilqg.XFSFull);
    XCOMFTorsoTDec = cat(3, XCOMFTorsoTDec, ilqg.XCOMFTorsoTDec);
    XCOMTTorsoFDec = cat(3, XCOMTTorsoFDec, ilqg.XCOMTTorsoFDec);
    
    UJoint = cat(3, UJoint, ilqg.UJoint);
    UCOMFS = cat(3, UCOMFS, ilqg.UCOMFS);
    UTorsoTS = cat(3, UTorsoTS, ilqg.UTorsoTS);
    UCOMTS = cat(3, UCOMTS, ilqg.UCOMTS);
    UTorsoFS = cat(3, UTorsoFS, ilqg.UTorsoFS);
    UTSFull = cat(3, UTSFull, ilqg.UTSFull);
    UFSFull = cat(3, UFSFull, ilqg.UFSFull);
    UCOMFTorsoTDec = cat(3, UCOMFTorsoTDec, ilqg.UCOMFTorsoTDec);
    UCOMTTorsoFDec = cat(3, UCOMTTorsoFDec, ilqg.UCOMTTorsoFDec);
    
%     CJoint  = cat(3, CJoint, ilqg.CJoint);
%     CTorsoTS = cat(3, CTorsoTS, ilqg.CTorsoTS);
%     CCOMFS = cat(3, CCOMFS, ilqg.CCOMFS);
%     CTorsoFS = cat(3, CTorsoFS, ilqg.CTorsoFS);
%     CCOMTS = cat(3, CCOMTS, ilqg.CCOMTS);
    
    XCOMFClose = cat(3, XCOMFClose, ilqg.XCOMFClose);
    XTorsoTClose = cat(3, XTorsoTClose, ilqg.XTorsoTClose);
    XCOMTClose = cat(3, XCOMTClose, ilqg.XCOMTClose);
    XTorsoFClose = cat(3, XTorsoFClose, ilqg.XTorsoFClose);
    XFFullClose = cat(3, XFFullClose, ilqg.XFFullClose);
    XTFullClose = cat(3, XTFullClose, ilqg.XTFullClose);
    
    UCOMFClose = cat(3, UCOMFClose, ilqg.UCOMFClose);
    UTorsoTClose = cat(3, UTorsoTClose, ilqg.UTorsoTClose);
    UCOMTClose = cat(3, UCOMTClose, ilqg.UCOMTClose);
    UTorsoFClose = cat(3, UTorsoFClose, ilqg.UTorsoFClose);
    UFFullClose = cat(3, UFFullClose, ilqg.UFFullClose);
    UTFullClose = cat(3, UTFullClose, ilqg.UTFullClose);
    
    KCOMFClose = cat(4, KCOMFClose, ilqg.KCOMFClose);
    KTorsoTClose = cat(4, KTorsoTClose, ilqg.KTorsoTClose);
    KCOMTClose = cat(4, KCOMTClose, ilqg.KCOMTClose);
    KTorsoFClose = cat(4, KTorsoFClose, ilqg.KTorsoFClose);
    KFFullClose = cat(4, KFFullClose, ilqg.KFFullClose);
    KTFullClose = cat(4, KTFullClose, ilqg.KTFullClose);
end

[x_starts, ia_x_starts, ~] = unique(x_starts', 'rows');
x_starts = x_starts';
% subset_starts = [1, 2, 5, 6, 7, 8, 11, 12, 19, 20, 23, 24, 25, 26, 29, 30];
% x_starts = x_starts(:, subset_starts);
% ia_x_starts = ia_x_starts(subset_starts);

XJoint = XJoint(:,:, ia_x_starts);
XCOMFS = XCOMFS(:,:, ia_x_starts);
XTorsoTS = XTorsoTS(:,:, ia_x_starts);
XCOMTS = XCOMTS(:,:, ia_x_starts);
XTorsoFS = XTorsoFS(:,:, ia_x_starts);
XTSFull = XTSFull(:,:, ia_x_starts);
XFSFull = XFSFull(:,:, ia_x_starts);
XCOMFTorsoTDec = XCOMFTorsoTDec(:,:, ia_x_starts);
XCOMTTorsoFDec = XCOMTTorsoFDec(:,:, ia_x_starts);

UJoint = UJoint(:,:, ia_x_starts);
UCOMFS = UCOMFS(:,:, ia_x_starts);
UTorsoTS = UTorsoTS(:,:, ia_x_starts);
UCOMTS = UCOMTS(:,:, ia_x_starts);
UTorsoFS = UTorsoFS(:,:, ia_x_starts);
UTSFull = UTSFull(:,:, ia_x_starts);
UFSFull = UFSFull(:,:, ia_x_starts);
UCOMFTorsoTDec = UCOMFTorsoTDec(:,:, ia_x_starts);
UCOMTTorsoFDec = UCOMTTorsoFDec(:,:, ia_x_starts);

% CJoint  = CJoint(:,:, ia_x_starts);
% CTorsoTS = CTorsoTS(:,:, ia_x_starts);
% CCOMFS = CCOMFS(:,:, ia_x_starts);
% CTorsoFS = CTorsoFS(:,:, ia_x_starts);
% CCOMTS = CCOMTS(:,:, ia_x_starts);

XCOMFClose = XCOMFClose(:,:, ia_x_starts);
XTorsoTClose = XTorsoTClose(:,:, ia_x_starts);
XCOMTClose = XCOMTClose(:,:, ia_x_starts);
XTorsoFClose = XTorsoFClose(:,:, ia_x_starts);
XFFullClose = XFFullClose(:,:, ia_x_starts);
XTFullClose = XTFullClose(:,:, ia_x_starts);

UCOMFClose = UCOMFClose(:,:, ia_x_starts);
UTorsoTClose = UTorsoTClose(:,:, ia_x_starts);
UCOMTClose = UCOMTClose(:,:, ia_x_starts);
UTorsoFClose = UTorsoFClose(:,:, ia_x_starts);
UFFullClose = UFFullClose(:,:, ia_x_starts);
UTFullClose = UTFullClose(:,:, ia_x_starts);

KCOMFClose = KCOMFClose(:,:,:, ia_x_starts);
KTorsoTClose = KTorsoTClose(:,:,:, ia_x_starts);
KCOMTClose = KCOMTClose(:,:,:, ia_x_starts);
KTorsoFClose = KTorsoFClose(:,:,:, ia_x_starts);
KFFullClose = KFFullClose(:,:,:, ia_x_starts);
KTFullClose = KTFullClose(:,:,:, ia_x_starts);

DP_dir = '../../Biped2D/data/';
DP_filename = 'BipedDecomposed.mat';
DP_joint_filename = 'Biped2D_joint.mat';
DP_extreme_filename = 'BipedDecomposedExtreme.mat';

load(strcat(DP_dir, DP_filename), 'V_T_Torso_second', 'F_COM_first_', 'T_Torso_second', ...
                                  'V_F_COM_second', 'T_Torso_first_', 'F_COM_second', ...
                                  'V_F_Torso_second', 'T_COM_first_', 'F_Torso_second', ...
                                  'V_T_COM_second', 'F_Torso_first_', 'T_COM_second', ...
                                  'V_F_COM_T_Torso', 'V_T_COM_F_Torso', ...
                                  'grid_l', 'grid_a', 'grid_x_dot', 'grid_z_dot', 'grid_th', 'grid_th_dot', ...
                                  'limits');
load(strcat(DP_dir, DP_joint_filename), 'V_joint', 'policy_joint');
load(strcat(DP_dir, DP_extreme_filename), 'paramsDec1S', 'paramsDec1F', 'paramsDec2S', 'paramsDec2F');

%% Estimate DDP values
J_DDP = nan(size(x_starts,2), 9);

DDP_finalerr = nan(size(x_starts,2), 9);
convergence_tol = 5e-2;

disp('DDP');
for jj=1:1:size(x_starts, 2)
    
    J_DDP(jj, 1)         = J_Biped2D(sys, XJoint(:,:, jj), UJoint(:,:, jj));
    lastState            = XJoint(:, end, jj);
    DDP_finalerr(jj, 1)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [3;4];
    sys.U_DIMS_FIXED     = [1;2];
    cost                 = l_Biped2DSecond(sys, XTorsoTS(:,:, jj), UTorsoTS(:,:, jj), ...
                                           UCOMFClose(:,:, jj), KCOMFClose(:,:,:, jj), XCOMFClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 2)         = sum(cost.*discount*sys.dt, 2);
    lastState            = XTorsoTS(:, end, jj);
    DDP_finalerr(jj, 2)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [1;2];
    sys.U_DIMS_FIXED     = [3;4];
    cost                 = l_Biped2DSecond(sys, XCOMFS(:,:, jj), UCOMFS(:,:, jj), ...
                                           UTorsoTClose(:,:, jj), KTorsoTClose(:,:,:, jj), XTorsoTClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 3)         = sum(cost.*discount*sys.dt, 2);
    lastState            = XCOMFS(:, end, jj);
    DDP_finalerr(jj, 3)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [1;2];
    sys.U_DIMS_FIXED     = [3;4];
    cost                 = l_Biped2DSecond(sys, XTorsoFS(:,:, jj), UTorsoFS(:,:, jj), ...
                                           UCOMTClose(:,:, jj), KCOMTClose(:,:,:, jj), XCOMTClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 4)         = sum(cost.*discount*sys.dt, 2);
    lastState            = XTorsoFS(:, end, jj);
    DDP_finalerr(jj, 4)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [3;4];
    sys.U_DIMS_FIXED     = [1;2];
    cost                 = l_Biped2DSecond(sys, XCOMTS(:,:, jj), UCOMTS(:,:, jj), ...
                                           UTorsoFClose(:,:, jj), KTorsoFClose(:,:,:, jj), XTorsoFClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 5)         = sum(cost.*discount*sys.dt, 2);
    lastState            = XCOMTS(:, end, jj);
    DDP_finalerr(jj, 5)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [3;4];
    sys.U_DIMS_FIXED     = [1;2];
    cost                 = l_Biped2DSecond(sys, XTSFull(:,:, jj), UTSFull(:,:, jj), ...
                                           UFFullClose(:,:, jj), KFFullClose(:,:,:, jj), XFFullClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 6)         = sum(cost.*discount*sys.dt, 2);
    lastState            = XTSFull(:, end, jj);
    DDP_finalerr(jj, 6)  = (norm(lastState - sys.goal));
    
    sys.U_DIMS_FREE      = [1;2];
    sys.U_DIMS_FIXED     = [3;4];
    cost                 = l_Biped2DSecond(sys, XFSFull(:,:, jj), UFSFull(:,:, jj), ...
                                           UTFullClose(:,:, jj), KTFullClose(:,:,:, jj), XTFullClose(:,:, jj));
    discount             = sys.gamma_.^linspace(0, length(cost)-1, length(cost));
    J_DDP(jj, 7)         = sum(cost.*discount*sys.dt, 2);
    lastState            = XFSFull(:, end, jj);
    DDP_finalerr(jj, 7)  = (norm(lastState - sys.goal));
    
    J_DDP(jj, 8)         = J_Biped2D(sys, XCOMFTorsoTDec(:,:, jj), UCOMFTorsoTDec(:,:, jj));
    lastState            = XCOMFTorsoTDec(:, end, jj);
    DDP_finalerr(jj, 8)  = (norm(lastState - sys.goal));
    
    J_DDP(jj, 9)         = J_Biped2D(sys, XCOMTTorsoFDec(:,:, jj), UCOMTTorsoFDec(:,:, jj));
    lastState            = XCOMTTorsoFDec(:, end, jj);
    DDP_finalerr(jj, 9)  = (norm(lastState - sys.goal));
    jj
end
DDP_converged = (DDP_finalerr < convergence_tol);

%% DP
disp('DP');
sys.U_DIMS_FREE = [1;2;3;4];
sys.U_DIMS_FIXED = [];
sys.X_DIMS_FREE = [1;2;3;4;5;6];
sys.X_DIMS_FIXED = [];

J_DP_rollout = nan(size(x_starts, 2), 9);
J_DP_longhorz_rollout = nan(size(x_starts, 2), 9);
DP_finalerr = nan(size(x_starts, 2), 9);
DP_longhorz_finalerr = nan(size(x_starts, 2), 9);
T = 4;
NUM_STEPS = round(T / sys.dt);
T_longhorz = 8;
NUM_STEPS_LONGHORZ = round(T_longhorz / sys.dt);

disp('Joint');
[trajXJoint, trajUJoint, J_DP_rollout(:,1)] = rolloutDPTrajOde(policy_joint, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXJointErr = trajXJoint(:, end, :) - sys.goal;
DP_finalerr(:, 1) = reshape(vecnorm(trajXJointErr), size(x_starts,2), 1);

[trajXJoint, trajUJoint, J_DP_longhorz_rollout(:,1)] = rolloutDPTrajOde(policy_joint, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXJointErr = trajXJoint(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 1) = reshape(vecnorm(trajXJointErr), size(x_starts,2), 1);

disp('Dec 1');
policyCOMFTorsoT(:,:,:,:,:,:,1:2) = F_COM_first_;
policyCOMFTorsoT(:,:,:,:,:,:,3:4) = T_Torso_second;

[trajXCFPT, trajUCFPT, J_DP_rollout(:,2)] = rolloutDPTrajOde(policyCOMFTorsoT, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXCFPTErr = trajXCFPT(:, end, :) - sys.goal;
DP_finalerr(:, 2) = reshape(vecnorm(trajXCFPTErr), size(x_starts,2), 1);

[trajXCFPT, trajUCFPT, J_DP_longhorz_rollout(:,2)] = rolloutDPTrajOde(policyCOMFTorsoT, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXCFPTErr = trajXCFPT(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 2) = reshape(vecnorm(trajXCFPTErr), size(x_starts,2), 1);

disp('Dec 2');
policyTorsoTCOMF(:,:,:,:,:,:,1:2) = F_COM_second;
policyTorsoTCOMF(:,:,:,:,:,:,3:4) = T_Torso_first_;

[trajXPTCF, trajUPTCF, J_DP_rollout(:,3)] = rolloutDPTrajOde(policyTorsoTCOMF, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXPTCFErr = trajXPTCF(:, end, :) - sys.goal;
DP_finalerr(:, 3) = reshape(vecnorm(trajXPTCFErr), size(x_starts,2), 1);

[trajXPTCF, trajUPTCF, J_DP_longhorz_rollout(:,3)] = rolloutDPTrajOde(policyTorsoTCOMF, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXPTCFErr = trajXPTCF(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 3) = reshape(vecnorm(trajXPTCFErr), size(x_starts,2), 1);

disp('Dec 3');
policyCOMTTorsoF(:,:,:,:,:,:,1:2) = F_Torso_second;
policyCOMTTorsoF(:,:,:,:,:,:,3:4) = T_COM_first_;

[trajXCTPF, trajUCTPF, J_DP_rollout(:,4)] = rolloutDPTrajOde(policyCOMTTorsoF, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXCTPFErr = trajXCTPF(:, end, :) - sys.goal;
DP_finalerr(:, 4) = reshape(vecnorm(trajXCTPFErr), size(x_starts,2), 1);

[trajXCTPF, trajUCTPF, J_DP_longhorz_rollout(:,4)] = rolloutDPTrajOde(policyCOMTTorsoF, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXCTPFErr = trajXCTPF(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 4) = reshape(vecnorm(trajXCTPFErr), size(x_starts,2), 1);

disp('Dec 4');
policyTorsoFCOMT(:,:,:,:,:,:,1:2) = F_Torso_first_;
policyTorsoFCOMT(:,:,:,:,:,:,3:4) = T_COM_second;
[trajXPFCT, trajUPFCT, J_DP_rollout(:,5)] = rolloutDPTrajOde(policyTorsoFCOMT, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXPFCTErr = trajXPFCT(:, end, :) - sys.goal;
DP_finalerr(:, 5) = reshape(vecnorm(trajXPFCTErr), size(x_starts,2), 1);

[trajXPFCT, trajUPFCT, J_DP_longhorz_rollout(:,5)] = rolloutDPTrajOde(policyTorsoFCOMT, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXPFCTErr = trajXPFCT(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 5) = reshape(vecnorm(trajXPFCTErr), size(x_starts,2), 1);

disp('Dec 5');
policyFullFFullT(:,:,:,:,:,:,1:2) = paramsDec2F.policy;
policyFullFFullT(:,:,:,:,:,:,3:4) = paramsDec2S.policy;
[trajXFFFT, trajUFFFT, J_DP_rollout(:,6)] = rolloutDPTrajOde(policyFullFFullT, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXFFFTErr = trajXFFFT(:, end, :) - sys.goal;
DP_finalerr(:, 6) = reshape(vecnorm(trajXFFFTErr), size(x_starts,2), 1);

[trajXFFFT, trajUFFFT, J_DP_longhorz_rollout(:,6)] = rolloutDPTrajOde(policyFullFFullT, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXFFFTErr = trajXFFFT(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 6) = reshape(vecnorm(trajXFFFTErr), size(x_starts,2), 1);

disp('Dec 6');
policyFullTFullF(:,:,:,:,:,:,1:2) = paramsDec1S.policy;
policyFullTFullF(:,:,:,:,:,:,3:4) = paramsDec1F.policy;
[trajXFTFF, trajUFTFF, J_DP_rollout(:,7)] = rolloutDPTrajOde(policyFullTFullF, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXFTFFErr = trajXFTFF(:, end, :) - sys.goal;
DP_finalerr(:, 7) = reshape(vecnorm(trajXFTFFErr), size(x_starts,2), 1);

[trajXFTFF, trajUFTFF, J_DP_longhorz_rollout(:,7)] = rolloutDPTrajOde(policyFullTFullF, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXFTFFErr = trajXFTFF(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 7) = reshape(vecnorm(trajXFTFFErr), size(x_starts,2), 1);

disp('Dec 7');
policyCOMFTorsoTDec(:,:,:,:,:,:,1:2) = F_COM_first_;
policyCOMFTorsoTDec(:,:,:,:,:,:,3:4) = T_Torso_first_;

[trajXCFPTD, trajUCFPTD, J_DP_rollout(:,8)] = rolloutDPTrajOde(policyCOMFTorsoTDec, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXCFPTDErr = trajXCFPTD(:, end, :) - sys.goal;
DP_finalerr(:, 8) = reshape(vecnorm(trajXCFPTDErr), size(x_starts,2), 1);

[trajXCFPTD, trajUCFPTD, J_DP_longhorz_rollout(:,8)] = rolloutDPTrajOde(policyCOMFTorsoTDec, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXCFPTDErr = trajXCFPTD(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 8) = reshape(vecnorm(trajXCFPTDErr), size(x_starts,2), 1);

disp('Dec 8');
policyCOMTTorsoFDec(:,:,:,:,:,:,1:2) = F_Torso_first_;
policyCOMTTorsoFDec(:,:,:,:,:,:,3:4) = T_COM_first_;

[trajXCTPFD, trajUCTPFD, J_DP_rollout(:,9)] = rolloutDPTrajOde(policyCOMTTorsoFDec, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T);
trajXCTPFDErr = trajXCTPFD(:, end, :) - sys.goal;
DP_finalerr(:, 9) = reshape(vecnorm(trajXCTPFDErr), size(x_starts,2), 1);

[trajXCTPFD, trajUCTPFD, J_DP_longhorz_rollout(:,9)] = rolloutDPTrajOde(policyCOMTTorsoFDec, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, T_longhorz);
trajXCTPFDErr = trajXCTPFD(:, end, :) - sys.goal;
DP_longhorz_finalerr(:, 9) = reshape(vecnorm(trajXCTPFDErr), size(x_starts,2), 1);

DP_converged = DP_finalerr < convergence_tol;
DP_longhorz_converged = DP_longhorz_finalerr < convergence_tol;

%% DP values for starts
disp('DP');
J_DP = nan(size(x_starts, 2), 9);
J_DP(:, 1) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, V_joint, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');
                 
J_DP(:, 2) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, V_t_pole_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');

J_DP(:, 3) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, V_f_cart_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');

J_DP(:, 4) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, V_f_pole_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');

J_DP(:, 5) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, V_t_cart_second, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');
                     
J_DP(:, 6) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, paramsDec2S.V, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');

J_DP(:, 7) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, paramsDec1S.V, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');

J_DP(:, 8) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, V_f_cart_t_pole, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');

J_DP(:, 9) = interpn(grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot, V_t_cart_f_pole, ...
                     x_starts(1,:)', x_starts(2,:)', x_starts(3,:)', x_starts(4,:)', x_starts(5,:)', x_starts(6,:)');

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

%% Compare trajectories
linew = 2;
% trajIDs = [1, 2, 5, 7, 11, 12, 13, 18, 19, 24, 25, 26, 31, 32, 33, 38, 39, 44, 45, 46, 49, 50, 51, 52, 55, 56];
trajIDs = find(~DP_converged(:,1));
% for jj=1:1:size(x_starts, 2)
for ii=1:1:length(trajIDs)
    jj = trajIDs(ii);
    
    figure;
    subplot(2,1,1);
    hold on;
    % plot DDP trajectories
    plot(XJoint(1,:,jj),XJoint(2,:,jj),'Color',[0.4660 0.6740 0.1880], 'LineWidth', linew);
    plot(XTorsoTS(1,:,jj),XTorsoTS(2,:,jj),'Color',[0 0.4470 0.7410], 'LineWidth', linew);
    plot(XCOMFS(1,:,jj),XCOMFS(2,:,jj),'Color',[0.8500 0.3250 0.0980], 'LineWidth', linew);
%     plot(XTorsoFS(1,:,jj),XTorsoFS(2,:,jj),'Color',[0.9290 0.6940 0.1250], 'LineWidth', linew);
%     plot(XCOMTS(1,:,jj),XCOMTS(2,:,jj),'Color',[0.4940 0.1840 0.5560], 'LineWidth', linew);
    plot(XCOMFTorsoTDec(1,:,jj), XCOMFTorsoTDec(2,:,jj), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', linew);
%     plot(XCOMTTorsoFDec(1,:,jj), XCOMTTorsoFDec(2,:,jj), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', linew);
    % plot DP trajectories
    plot(trajXJoint(1,:,jj), trajXJoint(2,:,jj),'Color',[0.4660 0.6740 0.1880], 'LineStyle', '-.', 'LineWidth', linew);
    plot(trajXCFPT(1,:,jj), trajXCFPT(2,:,jj),'Color',[0 0.4470 0.7410], 'LineStyle', '-.', 'LineWidth', linew);
    plot(trajXPTCF(1,:,jj), trajXPTCF(2,:,jj),'Color',[0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', linew);
%     plot(trajXCTPF(1,:,jj), trajXCTPF(2,:,jj),'Color',[0.9290 0.6940 0.1250], 'LineStyle', '-.', 'LineWidth', linew);
%     plot(trajXPFCT(1,:,jj), trajXPFCT(2,:,jj),'Color',[0.4940 0.1840 0.5560], 'LineStyle', '-.', 'LineWidth', linew);
    plot(trajXCFPTD(1,:,jj), trajXCFPTD(2,:,jj), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-.', 'LineWidth', linew);
%     plot(trajXCTPFD(1,:,jj), trajXCTPFD(2,:,jj), 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '-.', 'LineWidth', linew);
    
    legend(["Joint (DDP)", "F - COM, T - Both (DDP)", "T - Torso, F - Both (DDP)", ... % "T - COM, F - Both (DDP)", "F - Torso, T - Both (DDP)", ...
            "F - COM, T - Torso (DDP)", ... % "T - COM, F - Torso (DDP)", ...
            "Joint (DP)", "F - COM, T - Both (DP)", "T - Torso, F - Both (DP)", ... % "T - COM, F - Both (DP)", "F - Torso, T - Both (DP)", ...
            "F - COM, T - Torso (DP)", ... %"T - COM, F - Torso (DP)", ...
            ]);
    scatter(XJoint(1,1,jj),XJoint(2,1,jj), 20, [1,0,0]);
    xlabel('x');
    ylabel('x-dot');
    xlim([-1,1])
    ylim([-1.5,1.5])
    hold off;
    
    subplot(2,1,2);
    hold on;
    % plot DDP trajectories
    plot(XJoint(3,:,jj), XJoint(4,:,jj),'Color',[0.4660 0.6740 0.1880], 'LineWidth', linew);
    plot(XTorsoTS(3,:,jj), XTorsoTS(4,:,jj),'Color',[0 0.4470 0.7410], 'LineWidth', linew);
    plot(XCOMFS(3,:,jj), XCOMFS(4,:,jj),'Color',[0.8500 0.3250 0.0980], 'LineWidth', linew);
%     plot(XTorsoFS(3,:,jj), XTorsoFS(4,:,jj),'Color',[0.9290 0.6940 0.1250], 'LineWidth', linew);
%     plot(XCOMTS(3,:,jj), XCOMTS(4,:,jj),'Color',[0.4940 0.1840 0.5560], 'LineWidth', linew);
    plot(XCOMFTorsoTDec(3,:,jj), XCOMFTorsoTDec(4,:,jj), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', linew);
%     plot(XCOMTTorsoFDec(3,:,jj), XCOMTTorsoFDec(4,:,jj), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', linew);
    % plot DP trajectories
    plot(trajXJoint(3,:,jj), trajXJoint(4,:,jj),'Color',[0.4660 0.6740 0.1880], 'LineStyle', '-.', 'LineWidth', linew);
    plot(trajXCFPT(3,:,jj), trajXCFPT(4,:,jj),'Color',[0 0.4470 0.7410], 'LineStyle', '-.', 'LineWidth', linew);
    plot(trajXPTCF(3,:,jj), trajXPTCF(4,:,jj),'Color',[0.8500 0.3250 0.0980], 'LineStyle', '-.', 'LineWidth', linew);
%     plot(trajXCTPF(3,:,jj), trajXCTPF(4,:,jj),'Color',[0.9290 0.6940 0.1250], 'LineStyle', '-.', 'LineWidth', linew);
%     plot(trajXPFCT(3,:,jj), trajXPFCT(4,:,jj),'Color',[0.4940 0.1840 0.5560], 'LineStyle', '-.', 'LineWidth', linew);
    plot(trajXCFPTD(3,:,jj), trajXCFPTD(4,:,jj), 'Color', [0.3010 0.7450 0.9330], 'LineStyle', '-.', 'LineWidth', linew);
%     plot(trajXCTPFD(3,:,jj), trajXCTPFD(4,:,jj), 'Color', [0.6350 0.0780 0.1840], 'LineStyle', '-.', 'LineWidth', linew);
    
    legend(["Joint (DDP)", "F - COM, T - Both (DDP)", "T - Torso, F - Both (DDP)", ... % "T - COM, F - Both (DDP)", "F - Torso, T - Both (DDP)", ...
            "F - COM, T - Torso (DDP)", ... % "T - COM, F - Torso (DDP)", ...
            "Joint (DP)", "F - COM, T - Both (DP)", "T - Torso, F - Both (DP)", ... % "T - COM, F - Both (DP)", "F - Torso, T - Both (DP)", ...
            "F - COM, T - Torso (DP)", ... %"T - COM, F - Torso (DP)", ...
            ]);
    scatter(XJoint(3,1,jj),XJoint(4,1,jj),20,[1,0,0]);
    xlabel('theta');
    ylabel('theta-dot');
    xlim([pi, 3*pi])
    ylim([-3, 3])
    hold off;
end

%% DP Value function comparison within state bounds

state_bounds = [-0.92, 1;
                pi/2 + 0.3, pi/2 + 0.4;
                -0.1, 0.1;
                -0.4, -0.3;
                -0.4, 0.4;
                -0.2, 0.2];

valid_range = ((grid_l >= state_bounds(1,1)) & (grid_l <= state_bounds(1,2)) ...
                & (grid_a >= state_bounds(2,1)) & (grid_a <= state_bounds(2,2)) ...
                & (grid_x_dot >= state_bounds(3,1)) & (grid_x_dot <= state_bounds(3,2)) ...
                & (grid_z_dot >= state_bounds(4,1)) & (grid_z_dot <= state_bounds(4,2)) ...
                & (grid_th >= state_bounds(5,1)) & (grid_th <= state_bounds(5,2)) ...
                & (grid_th_dot >= state_bounds(6,1)) & (grid_th_dot <= state_bounds(6,2)));

DP_avg_err = [sum(abs(V_joint(valid_range) - V_T_Torso_second(valid_range)), 'all'), ... 
              sum(abs(V_joint(valid_range) - V_F_COM_second(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_F_Torso_second(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_T_COM_second(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - paramsDec2S.V(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - paramsDec1S.V(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_F_COM_T_Torso(valid_range)), 'all'), ...
              sum(abs(V_joint(valid_range) - V_T_COM_F_Torso(valid_range)), 'all')];

[~, DP_avg_order] = sort(DP_avg_err);

%% Inidividual test

% [tXJ, tUJ, J, JW] = rolloutDPTrajOde(policy_joint, x_starts(:,1), sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T);
% [tXJL, tUJL, JL, JWL] = rolloutDPTrajOde(policy_joint, x_starts(:,1), sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, T_longhorz);
% [tXJLL, tUJLL, JLL, JWLL] = rolloutDPTrajOde(policy_joint, x_starts(:,1), sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, 2*T_longhorz);
% [tXJLLL, tUJLLL, JLLL, JWLLL] = rolloutDPTrajOde(policy_joint, x_starts(:,1), sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, 4*T_longhorz);
% [tXJLLLL, tUJLLLL, JLLLL, JWLLLL] = rolloutDPTrajOde(policy_joint, x_starts(:,1), sys, limits, grid_x, grid_x_dot, grid_xP, grid_xP_dot, 32*T_longhorz);

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
        J_WrapAround(jj, 1) = J_Biped2D_WrapAround(sys, trajX(:,:, jj), trajU(:,:, jj));
        J(jj, 1) = J_Biped2D(sys, trajX(:,:, jj), trajU(:,:, jj));
    end
end

function [trajX, trajU, J] = rolloutDPTrajOde(policy, x_starts, sys, limits, grid_l, grid_a, grid_x_dot, grid_xP_dot, grid_th, grid_th_dot, T)
    
    NUM_CTRL = round(T / sys.dt);
    tspan = linspace(0, T, NUM_CTRL+1);
    opts = odeset('AbsTol', 1e-4, 'RelTol', 1e-4);
    trajX = nan(6, NUM_CTRL+1, size(x_starts, 2));
    trajU = nan(4, NUM_CTRL, size(x_starts, 2));
    J = zeros(size(x_starts, 2), 1);
    P1 = griddedInterpolant(grid_l, grid_a, grid_x_dot, grid_xP_dot, grid_th, grid_th_dot, policy(:,:,:,:,:,:,1));
    P2 = griddedInterpolant(grid_l, grid_a, grid_x_dot, grid_xP_dot, grid_th, grid_th_dot, policy(:,:,:,:,:,:,2));
    P3 = griddedInterpolant(grid_l, grid_a, grid_x_dot, grid_xP_dot, grid_th, grid_th_dot, policy(:,:,:,:,:,:,3));
    P4 = griddedInterpolant(grid_l, grid_a, grid_x_dot, grid_xP_dot, grid_th, grid_th_dot, policy(:,:,:,:,:,:,4));
    for jj = 1:1:size(x_starts, 2)
        [~, X] = ode113(@(t,y) biped2d_dyn_gridbased(t, y, P1, P2, P3, P4, limits, sys), tspan, x_starts(:, jj), opts);
        trajX(:,:, jj) = X';
        X_ = [min(limits(1, 2), max(limits(1, 1), trajX(1, :, jj)));
              min(limits(2, 2), max(limits(2, 1), trajX(2, :, jj)));
              min(limits(3, 2), max(limits(3, 1), trajX(3, :, jj)));
              min(limits(4, 2), max(limits(4, 1), trajX(4, :, jj)));
              min(limits(5, 2), max(limits(5, 1), trajX(5, :, jj)));
              min(limits(6, 2), max(limits(6, 1), trajX(6, :, jj)))];
        
        trajU(1, :, jj) = P1(X_(1,1:(end-1)), X_(2,1:(end-1)), X_(3,1:(end-1)), X_(4,1:(end-1)), X_(5,1:(end-1)), X_(6,1:(end-1)));
        trajU(2, :, jj) = P2(X_(1,1:(end-1)), X_(2,1:(end-1)), X_(3,1:(end-1)), X_(4,1:(end-1)), X_(5,1:(end-1)), X_(6,1:(end-1)));
        trajU(3, :, jj) = P3(X_(1,1:(end-1)), X_(2,1:(end-1)), X_(3,1:(end-1)), X_(4,1:(end-1)), X_(5,1:(end-1)), X_(6,1:(end-1)));
        trajU(4, :, jj) = P4(X_(1,1:(end-1)), X_(2,1:(end-1)), X_(3,1:(end-1)), X_(4,1:(end-1)), X_(5,1:(end-1)), X_(6,1:(end-1)));
        J(jj, 1) = J_Biped2D(sys, trajX(:,:, jj), trajU(:,:, jj));
        disp(strcat('start : ', num2str(jj)));
    end
    
end