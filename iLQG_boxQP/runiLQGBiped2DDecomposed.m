clear;
clc;
close all;
addpath('biped2d');

%% Common parameters
sys.full_DDP = false;   
sys.m = 72;
sys.I = 3;
sys.l0 = 1.05;
sys.g = 9.81;
sys.d = 0.2;
sys.df = 0.5;
sys.T = 5;
sys.dt = 0.001;
NUM_CTRL = round(sys.T / sys.dt);
sys.gamma_ = 0.997;
sys.Q = diag([100, 200, 2, 2, 1000, 10]);
sys.R = 0.001*eye(4);

lg = 0.96;
alpha1g = pi/2 + asin(sys.df/lg);
sys.goal = [lg; alpha1g; 0; 0; 0; 0];
sys.l_point = sys.goal;
sys.lims = [0, 1.5*sys.m*sys.g;
            0, 1.5*sys.m*sys.g;
            -0.125*sys.m*sys.g, 0.125*sys.m*sys.g;
            -0.125*sys.m*sys.g, 0.125*sys.m*sys.g];
% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 200;
Op.gamma_ = sys.gamma_;
Op.Alpha = [1];

% Define starts
com_pos = [0.85, 0.85, 0.9, 0.9, 0.95, 0.95;
            0.1,  0.4, 0.1, 0.4,  0.1, 0.4];
com_pos(2,:) = pi/2 + com_pos(2,:);
com_vel = [ 0.1, 0.05, 0.02,  0.1, 0.05, 0.02;
           -0.3, -0.3, -0.3, -0.4, -0.4, -0.4];
theta_starts = [-0.1,  -0.05, 0.05, 0.1;
                   0,      0,    0,   0];

x_starts = nan(6, size(com_pos,2)*size(com_vel,2)*size(theta_starts,2));
ind = 1;
for ii=1:1:size(com_pos, 2)
    for jj=1:1:size(com_vel, 2)
        for kk=1:1:size(theta_starts, 2)
            start = [com_pos(:, ii); com_vel(:, jj); theta_starts(:, kk)];
            x_starts(:, ind) = start;
            ind = ind + 1;
        end
    end
end

x_starts = [0.9; pi/2 + 0.3; -0.2; -0.4; 0.3; 0];

%% COM - F, Torso - T
% COM First
sys_COMFF = sys;
sys_COMFF.U_DIMS_FREE = [1; 2];
sys_COMFF.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMFF.U_DIMS_FIXED(sys_COMFF.U_DIMS_FREE) = [];
sys_COMFF.X_DIMS_FREE = [1;2;3;4];
sys_COMFF.X_DIMS_FIXED = linspace(1,6,6)';
sys_COMFF.X_DIMS_FIXED(sys_COMFF.X_DIMS_FREE) = [];

Op.lims = sys_COMFF.lims(sys_COMFF.U_DIMS_FREE, :);
XCOMFF = zeros(length(sys_COMFF.X_DIMS_FREE) + length(sys_COMFF.X_DIMS_FIXED), NUM_CTRL+1, size(x_starts, 2));
UCOMFF = zeros(length(sys_COMFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMFF = zeros(length(sys_COMFF.U_DIMS_FREE), length(sys_COMFF.X_DIMS_FREE) + length(sys_COMFF.X_DIMS_FIXED), ...
               NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_COMFF = @(x, u, i) biped2d_dyn_first_cst(sys_COMFF, x, u, sys_COMFF.full_DDP);
    Uinit = zeros(length(sys_COMFF.U_DIMS_FREE), NUM_CTRL);
    [XCOMFF(sys_COMFF.X_DIMS_FREE,:, jj), ...
     UCOMFF(:,1:NUM_CTRL, jj), ...
     KCOMFF(:,sys_COMFF.X_DIMS_FREE,1:NUM_CTRL, jj)] = iLQG(dyn_COMFF, ...
                                                    x_starts(sys_COMFF.X_DIMS_FREE, jj), ...
                                                    Uinit, Op);
    XCOMFF(sys_COMFF.X_DIMS_FIXED,:, jj) = sys_COMFF.l_point(sys_COMFF.X_DIMS_FIXED) ...
                                              * ones(1, size(XCOMFF, 2));
end

% Torso second
sys_TorsoTS = sys;
sys_TorsoTS.U_DIMS_FREE = [3; 4];
sys_TorsoTS.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoTS.U_DIMS_FIXED(sys_TorsoTS.U_DIMS_FREE) = [];
sys_TorsoTS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_TorsoTS.X_DIMS_FIXED = [];

Op.lims = sys_TorsoTS.lims(sys_TorsoTS.U_DIMS_FREE, :);
XTorsoTS = zeros(length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoTS = zeros(length(sys_TorsoTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoTS = zeros(length(sys_TorsoTS.U_DIMS_FREE), length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XCOMFClose = nan(length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMFClose = zeros(length(sys_COMFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMFClose = zeros(length(sys_COMFF.U_DIMS_FREE), length(sys_TorsoTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_TorsoTS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_TorsoTS, x, u, k, K, xn, sys_TorsoTS.full_DDP);
    Uinit = zeros(length(sys_TorsoTS.U_DIMS_FREE), NUM_CTRL);
    [XTorsoTS(:,:, jj), ...
     UTorsoTS(:,1:NUM_CTRL, jj), ...
     KTorsoTS(:,:,1:NUM_CTRL, jj), ...
     UCOMFClose(:,:, jj), ...
     KCOMFClose(:,:,:, jj), ...
     XCOMFClose(:,:, jj)] = iLQGSecond(dyn_TorsoTS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UCOMFF(:,:, jj), ...
                                        KCOMFF(:,:,:, jj), ...
                                        XCOMFF(:,:, jj), ...
                                        Op);
end

%% Torso - T, COM - F
% Torso first
sys_TorsoTF = sys;
sys_TorsoTF.U_DIMS_FREE = [3; 4];
sys_TorsoTF.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoTF.U_DIMS_FIXED(sys_TorsoTF.U_DIMS_FREE) = [];
sys_TorsoTF.X_DIMS_FREE = [5; 6];
sys_TorsoTF.X_DIMS_FIXED = linspace(1,6,6)';
sys_TorsoTF.X_DIMS_FIXED(sys_TorsoTF.X_DIMS_FREE) = [];

Op.lims = sys_TorsoTF.lims(sys_TorsoTF.U_DIMS_FREE, :);
XTorsoTF = zeros(length(sys_TorsoTF.X_DIMS_FREE) + length(sys_TorsoTF.X_DIMS_FIXED), NUM_CTRL+1, size(x_starts, 2));
UTorsoTF = zeros(length(sys_TorsoTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoTF = zeros(length(sys_TorsoTF.U_DIMS_FREE), length(sys_TorsoTF.X_DIMS_FREE) + length(sys_TorsoTF.X_DIMS_FIXED), ...
                 NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_TorsoTF = @(x, u, i) biped2d_dyn_first_cst(sys_TorsoTF, x, u, sys_TorsoTF.full_DDP);
    Uinit = zeros(length(sys_TorsoTF.U_DIMS_FREE), NUM_CTRL);
    [XTorsoTF(sys_TorsoTF.X_DIMS_FREE,:, jj), ...
     UTorsoTF(:,1:NUM_CTRL, jj), ...
     KTorsoTF(:,sys_TorsoTF.X_DIMS_FREE,1:NUM_CTRL, jj)] = iLQG(dyn_TorsoTF, ...
                                                    x_starts(sys_TorsoTF.X_DIMS_FREE, jj), ...
                                                    Uinit, Op);
    XTorsoTF(sys_TorsoTF.X_DIMS_FIXED,:, jj) = sys_TorsoTF.l_point(sys_TorsoTF.X_DIMS_FIXED) ...
                                               * ones(1, size(XTorsoTF, 2));
end

% COM second
sys_COMFS = sys;
sys_COMFS.U_DIMS_FREE = [1; 2];
sys_COMFS.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMFS.U_DIMS_FIXED(sys_COMFS.U_DIMS_FREE) = [];
sys_COMFS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMFS.X_DIMS_FIXED = [];

Op.lims = sys_COMFS.lims(sys_COMFS.U_DIMS_FREE, :);
XCOMFS = zeros(length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMFS = zeros(length(sys_COMFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMFS = zeros(length(sys_COMFS.U_DIMS_FREE), length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTorsoTClose = nan(length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoTClose = zeros(length(sys_TorsoTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoTClose = zeros(length(sys_TorsoTF.U_DIMS_FREE), length(sys_COMFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_COMFS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_COMFS, x, u, k, K, xn, sys_COMFS.full_DDP);
    Uinit = zeros(length(sys_COMFS.U_DIMS_FREE), NUM_CTRL);
    [XCOMFS(:,:, jj), ...
     UCOMFS(:,1:NUM_CTRL, jj), ...
     KCOMFS(:,:,1:NUM_CTRL, jj), ...
     UTorsoTClose(:,:, jj), ...
     KTorsoTClose(:,:,:, jj), ...
     XTorsoTClose(:,:, jj)] = iLQGSecond(dyn_COMFS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UTorsoTF(:,:, jj), ...
                                        KTorsoTF(:,:,:, jj), ...
                                        XTorsoTF(:,:, jj), ...
                                        Op);
end

%% COM - T, Torso - F
% COM first
sys_COMTF = sys;
sys_COMTF.U_DIMS_FREE = [3; 4];
sys_COMTF.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMTF.U_DIMS_FIXED(sys_COMTF.U_DIMS_FREE) = [];
sys_COMTF.X_DIMS_FREE = [1;2;3;4];
sys_COMTF.X_DIMS_FIXED = linspace(1,6,6)';
sys_COMTF.X_DIMS_FIXED(sys_COMFF.X_DIMS_FREE) = [];

Op.lims = sys_COMTF.lims(sys_COMTF.U_DIMS_FREE, :);
XCOMTF = zeros(length(sys_COMTF.X_DIMS_FREE) + length(sys_COMTF.X_DIMS_FIXED), NUM_CTRL+1, size(x_starts, 2));
UCOMTF = zeros(length(sys_COMTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMTF = zeros(length(sys_COMTF.U_DIMS_FREE), length(sys_COMTF.X_DIMS_FREE) + length(sys_COMTF.X_DIMS_FIXED), ...
                NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_COMTF = @(x, u, i) biped2d_dyn_first_cst(sys_COMTF, x, u, sys_COMTF.full_DDP);
    Uinit = zeros(length(sys_COMTF.U_DIMS_FREE), NUM_CTRL);
    [XCOMTF(sys_COMTF.X_DIMS_FREE,:, jj), ...
     UCOMTF(:,1:NUM_CTRL, jj), ...
     KCOMTF(:,sys_COMTF.X_DIMS_FREE,1:NUM_CTRL, jj)] = iLQG(dyn_COMTF, ...
                                                    x_starts(sys_COMTF.X_DIMS_FREE, jj), ...
                                                    Uinit, Op);
    XCOMTF(sys_COMTF.X_DIMS_FIXED,:, jj) = sys_COMTF.l_point(sys_COMTF.X_DIMS_FIXED)...
                                              * ones(1, size(XCOMTF, 2));
end

% Torso second
sys_TorsoFS = sys;
sys_TorsoFS.U_DIMS_FREE = [1; 2];
sys_TorsoFS.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoFS.U_DIMS_FIXED(sys_TorsoFS.U_DIMS_FREE) = [];
sys_TorsoFS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_TorsoFS.X_DIMS_FIXED = [];

Op.lims = sys_TorsoFS.lims(sys_TorsoFS.U_DIMS_FREE, :);
XTorsoFS = zeros(length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoFS = zeros(length(sys_TorsoFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoFS = zeros(length(sys_TorsoFS.U_DIMS_FREE), length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XCOMTClose = nan(length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMTClose = zeros(length(sys_COMTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMTClose = zeros(length(sys_COMTF.U_DIMS_FREE), length(sys_TorsoFS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_TorsoFS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_TorsoFS, x, u, k, K, xn, sys_TorsoFS.full_DDP);
    Uinit = zeros(length(sys_TorsoFS.U_DIMS_FREE), NUM_CTRL);
    [XTorsoFS(:,:, jj), ...
     UTorsoFS(:,1:NUM_CTRL, jj), ...
     KTorsoFS(:,:,1:NUM_CTRL, jj), ...
     UCOMTClose(:,:, jj), ...
     KCOMTClose(:,:,:, jj), ...
     XCOMTClose(:,:, jj)] = iLQGSecond(dyn_TorsoFS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UCOMTF(:,:, jj), ...
                                        KCOMTF(:,:,:, jj), ...
                                        XCOMTF(:,:, jj), ...
                                        Op);
end

%% Torso - F, COM - T
% Torso first
sys_TorsoFF = sys;
sys_TorsoFF.U_DIMS_FREE = [1; 2];
sys_TorsoFF.U_DIMS_FIXED = linspace(1,4,4)';
sys_TorsoFF.U_DIMS_FIXED(sys_TorsoFF.U_DIMS_FREE) = [];
sys_TorsoFF.X_DIMS_FREE = [5; 6];
sys_TorsoFF.X_DIMS_FIXED = linspace(1,6,6)';
sys_TorsoFF.X_DIMS_FIXED(sys_TorsoFF.X_DIMS_FREE) = [];

Op.lims = sys_TorsoFF.lims(sys_TorsoFF.U_DIMS_FREE, :);
XTorsoFF = zeros(length(sys_TorsoFF.X_DIMS_FREE) + length(sys_TorsoFF.X_DIMS_FIXED), NUM_CTRL+1, size(x_starts, 2));
UTorsoFF = zeros(length(sys_TorsoFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoFF = zeros(length(sys_TorsoFF.U_DIMS_FREE), length(sys_TorsoFF.X_DIMS_FREE) + length(sys_TorsoFF.X_DIMS_FIXED), ...
                NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_TorsoFF = @(x, u, i) biped2d_dyn_first_cst(sys_TorsoFF, x, u, sys_TorsoFF.full_DDP);
    Uinit = zeros(length(sys_TorsoFF.U_DIMS_FREE), NUM_CTRL);
    [XTorsoFF(sys_TorsoFF.X_DIMS_FREE,:, jj), ...
     UTorsoFF(:,1:NUM_CTRL, jj), ...
     KTorsoFF(:,sys_TorsoFF.X_DIMS_FREE,1:NUM_CTRL, jj)] = iLQG(dyn_TorsoFF, ...
                                                    x_starts(sys_TorsoFF.X_DIMS_FREE, jj), ...
                                                    Uinit, Op);
    XTorsoFF(sys_TorsoFF.X_DIMS_FIXED,:, jj) = sys_TorsoFF.l_point(sys_TorsoFF.X_DIMS_FIXED)...
                                              * ones(1, size(XTorsoFF, 2));
end

% COM second
sys_COMTS = sys;
sys_COMTS.U_DIMS_FREE = [3; 4];
sys_COMTS.U_DIMS_FIXED = linspace(1,4,4)';
sys_COMTS.U_DIMS_FIXED(sys_COMTS.U_DIMS_FREE) = [];
sys_COMTS.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMTS.X_DIMS_FIXED = [];

Op.lims = sys_COMTS.lims(sys_COMTS.U_DIMS_FREE, :);
XCOMTS = zeros(length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMTS = zeros(length(sys_COMTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCOMTS = zeros(length(sys_COMTS.U_DIMS_FREE), length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
XTorsoFClose = nan(length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UTorsoFClose = zeros(length(sys_TorsoFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KTorsoFClose = zeros(length(sys_TorsoFF.U_DIMS_FREE), length(sys_COMTS.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_COMTS = @(x, u, k, K, xn, i) ...
               biped2d_dyn_second_cst(sys_COMTS, x, u, k, K, xn, sys_COMTS.full_DDP);
    Uinit = zeros(length(sys_COMTS.U_DIMS_FREE), NUM_CTRL);
    [XCOMTS(:,:, jj), ...
     UCOMTS(:,1:NUM_CTRL, jj), ...
     KCOMTS(:,:,1:NUM_CTRL, jj), ...
     UTorsoFClose(:,:, jj), ...
     KTorsoFClose(:,:,:, jj), ...
     XTorsoFClose(:,:, jj)] = iLQGSecond(dyn_COMTS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UTorsoFF(:,:, jj), ...
                                        KTorsoFF(:,:,:, jj), ...
                                        XTorsoFF(:,:, jj), ...
                                        Op);
end

%% Decoupled
% COM - F, Torso - T
sys_COMFTorsoTDec = sys;
sys_COMFTorsoTDec.U_DIMS_FREE = [1;2;3;4];
sys_COMFTorsoTDec.U_DIMS_FIXED = [];
sys_COMFTorsoTDec.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMFTorsoTDec.X_DIMS_FIXED = [];
XCOMFTorsoTDec = zeros(length(sys_COMFTorsoTDec.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMFTorsoTDec = zeros(length(sys_COMFTorsoTDec.U_DIMS_FREE), NUM_CTRL, size(x_starts, 2));
U_DIMS = cell(2,1);
U_DIMS{1,1} = [1;2];
U_DIMS{2,1} = [3;4];
for jj=1:1:size(x_starts, 2)
    dyn_COMFTorsoTDec = @(x, u) ...
               biped2d_dyn_first_cst(sys_COMFTorsoTDec, x, u, sys_COMFTorsoTDec.full_DDP);
    [XCOMFTorsoTDec(:,:, jj), UCOMFTorsoTDec(:,:, jj)] = ForwardPassDec([UCOMFF(:,:, jj); UTorsoTF(:,:, jj)], ...
                                                                        [KCOMFF(:,:,:, jj); KTorsoTF(:,:,:, jj)], ...
                                                                        [XCOMFF(:,:, jj); XTorsoTF(:,:, jj)], ...
                                                                        U_DIMS, ...
                                                                         x_starts(:, jj), dyn_COMFTorsoTDec);
    jj
end

% COM - T, Torso - F
sys_COMTTorsoFDec = sys;
sys_COMTTorsoFDec.U_DIMS_FREE = [1;2;3;4];
sys_COMTTorsoFDec.U_DIMS_FIXED = [];
sys_COMTTorsoFDec.X_DIMS_FREE = [1;2;3;4;5;6];
sys_COMTTorsoFDec.X_DIMS_FIXED = [];
XCOMTTorsoFDec = zeros(length(sys_COMTTorsoFDec.X_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
UCOMTTorsoFDec = zeros(length(sys_COMTTorsoFDec.U_DIMS_FREE), NUM_CTRL, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_COMTTorsoFDec = @(x, u) ...
               biped2d_dyn_first_cst(sys_COMTTorsoFDec, x, u, sys_COMTTorsoFDec.full_DDP);
    [XCOMTTorsoFDec(:,:, jj), UCOMTTorsoFDec(:,:, jj)] = ForwardPassDec([UTorsoFF(:,:, jj); UCOMTF(:,:, jj)], ...
                                                                        [KTorsoFF(:,:,:, jj); KCOMTF(:,:,:, jj)], ...
                                                                        [XTorsoFF(:,:, jj); XCOMTF(:,:, jj)], ...
                                                                        U_DIMS, ...
                                                                        x_starts(:, jj), dyn_COMTTorsoFDec);
    jj
end

%% Joint
tic;
sys_joint = sys;
sys_joint.U_DIMS_FREE = [1;2;3;4];
sys_joint.U_DIMS_FIXED = [];
sys_joint.X_DIMS_FREE = [1;2;3;4;5;6];
sys_joint.X_DIMS_FIXED = [];

Ops.lims = sys_joint.lims;
XJoint = zeros(4, NUM_CTRL+1, size(x_starts, 2));
UJoint = zeros(length(sys_joint.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KJoint = zeros(length(sys_joint.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));

for jj=1:1:size(x_starts, 2)
    dyn_joint = @(x, u, i) biped2d_dyn_first_cst(sys_joint, x, u, sys_joint.full_DDP);
    Uinit = zeros(length(sys_joint.U_DIMS_FREE), NUM_CTRL);
    [XJoint(sys_joint.X_DIMS_FREE,:, jj), ...
     UJoint(:,1:NUM_CTRL, jj), ...
     KJoint(:,sys_joint.X_DIMS_FREE,1:NUM_CTRL, jj)] = iLQG(dyn_joint, ...
                                                    x_starts(sys_joint.X_DIMS_FREE, jj), ...
                                                    Uinit, Op);
    XJoint(sys_joint.X_DIMS_FIXED,:, jj) = repmat(sys_joint.l_point(sys_joint.X_DIMS_FIXED), [1, size(XJoint,2)]);
    jj
end

status = "zeroinit_finallpoint_nowraparound";