clear;
clc;
close all;
addpath('cartpole');

%% Common parameters
sys.full_DDP = false;   
sys.mc = 5;
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81;
sys.T = 5;
sys.dt = 0.001;
NUM_CTRL = round(sys.T / sys.dt);
sys.gamma_ = 0.997;
sys.Q = [25, 0, 0, 0;
         0, 0.02, 0, 0;
         0, 0, 25, 0;
         0, 0, 0, 0.02];
sys.R = [0.001, 0;
         0, 0.001];
sys.goal = [0; 0; pi; 0];
sys.l_point = [0; 0; pi; 0];
sys.lims = [-9, 9;
            -9, 9];
% Optimization parameters
Op.lims  = sys.lims;
Op.maxIter = 200;
Op.gamma_ = sys.gamma_;
Op.Alpha = [0.25];

% Define starts
%cart_starts = [-1, -0.5, 0, 0.5, 1;
%                0,  0,   0, 0,   0];
%cart_starts = [-0.75, -0.5, 0, 0.5, 0.75;
%                0,     0,   0, 0,   0];
%cart_starts = [-0.4, -0.2, 0, 0.2, 0.4;
%                0,     0,  0, 0,   0];
cart_starts = [-0.4, -0.2, 0.2, 0.4;
                0,     0,  0,   0];

%pole_starts = [7*pi/4, 5*pi/4, 3*pi/4, pi/4, 0;
%               0,      0,      0,      0,    0];
%pole_starts = [pi/2, 2*pi/3, 5*pi/6, 7*pi/6, 4*pi/3, 3*pi/2;
%               0,      0,      0,      0,    0,      0];
%pole_starts = [2*pi/3, 3*pi/4, 5*pi/6, 7*pi/6, 5*pi/4, 4*pi/3;
%               0,      0,      0,      0,    0,      0];
pole_starts = [2*pi/3, 3*pi/4, 5*pi/4, 4*pi/3;
               0,      0,        0,        0];

x_starts = nan(4, size(cart_starts,2)*size(pole_starts,2));
for ii=1:1:size(cart_starts, 2)
    for jj=1:1:size(pole_starts, 2)
        x_starts(:, (ii-1)*size(pole_starts, 2) + jj) = [cart_starts(:, ii); pole_starts(:, jj)];
    end
end

%load('data/x_starts_notconverged.mat');
%x_starts = x_starts_notconverged;

%% Joint
disp('**** Joint ****');

sys_joint = sys;
sys_joint.U_DIMS_FREE = [1;2];
sys_joint.U_DIMS_FIXED = [];
sys_joint.X_DIMS_FREE = [1;2;3;4];
sys_joint.X_DIMS_FIXED = [];

Op.lims = sys_joint.lims;
XJoint = zeros(4, NUM_CTRL+1, size(x_starts, 2));
CJoint = zeros(1, NUM_CTRL+1, size(x_starts, 2));
UJoint = zeros(length(sys_joint.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KJoint = zeros(length(sys_joint.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
TraceJoint = cell(size(x_starts, 2), 1);
timeJoint = zeros(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    dyn_joint = @(x, u, i) cartpole_dyn_first_cst(sys_joint, x, u, sys_joint.full_DDP);
    Uinit = zeros(length(sys_joint.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XJoint(sys_joint.X_DIMS_FREE,:, jj), ...
     UJoint(:,1:NUM_CTRL, jj), ...
     KJoint(:,sys_joint.X_DIMS_FREE,1:NUM_CTRL, jj), ...
     ~, ~, CJoint(:,:, jj), TraceJoint{jj,1}] = iLQG(dyn_joint, ...
                                   x_starts(sys_joint.X_DIMS_FREE, jj), ...
                                   Uinit, Op);
    timeJoint(jj) = toc;
    XJoint(sys_joint.X_DIMS_FIXED,:, jj) = repmat(sys_joint.l_point(sys_joint.X_DIMS_FIXED), [1, size(XJoint,2)]);
    jj
end

%% Cart - F, Pole - T
disp('**** F - Cart, T - Both ****');
% Cart First
disp('F - Cart')
sys_CartFF = sys;
sys_CartFF.U_DIMS_FREE = [1];
sys_CartFF.U_DIMS_FIXED = linspace(1,2,2)';
sys_CartFF.U_DIMS_FIXED(sys_CartFF.U_DIMS_FREE) = [];
sys_CartFF.X_DIMS_FREE = [1;2];
sys_CartFF.X_DIMS_FIXED = linspace(1,4,4)';
sys_CartFF.X_DIMS_FIXED(sys_CartFF.X_DIMS_FREE) = [];

Op.lims = sys_CartFF.lims(sys_CartFF.U_DIMS_FREE, :);
XCartFF = zeros(4, NUM_CTRL+1, size(cart_starts, 2));
CCartFF = zeros(1, NUM_CTRL+1, size(cart_starts, 2));
UCartFF = zeros(length(sys_CartFF.U_DIMS_FREE), NUM_CTRL+1, size(cart_starts, 2));
KCartFF = zeros(length(sys_CartFF.U_DIMS_FREE), 4, NUM_CTRL+1, size(cart_starts, 2));
TraceCartFF = cell(size(cart_starts, 2), 1);
timeCartFF = zeros(size(cart_starts, 2), 1);
for jj=1:1:size(cart_starts, 2)
    dyn_CartFF = @(x, u, i) cartpole_dyn_first_cst(sys_CartFF, x, u, sys_CartFF.full_DDP);
    Uinit = zeros(length(sys_CartFF.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XCartFF(sys_CartFF.X_DIMS_FREE,:, jj), ...
     UCartFF(:,1:NUM_CTRL, jj), ...
     KCartFF(:,sys_CartFF.X_DIMS_FREE,1:NUM_CTRL, jj), ...
     ~, ~, CCartFF(:,:, jj), TraceCartFF{jj,1}] = iLQG(dyn_CartFF, ...
                                    x_starts(sys_CartFF.X_DIMS_FREE, jj), ...
                                    Uinit, Op);
    timeCartFF(jj) = toc;
    XCartFF(sys_CartFF.X_DIMS_FIXED,:, jj) = sys_CartFF.l_point(sys_CartFF.X_DIMS_FIXED)...
                                              * ones(1, size(XCartFF, 2));
    jj
end

% Pole second
disp('T - Both')
sys_PoleTS = sys;
sys_PoleTS.U_DIMS_FREE = [2];
sys_PoleTS.U_DIMS_FIXED = linspace(1,2,2)';
sys_PoleTS.U_DIMS_FIXED(sys_PoleTS.U_DIMS_FREE) = [];
sys_PoleTS.X_DIMS_FREE = [1;2;3;4];
sys_PoleTS.X_DIMS_FIXED = [];

Op.lims = sys_PoleTS.lims(sys_PoleTS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [1;2];
XPoleTS = zeros(4, NUM_CTRL+1, size(x_starts, 2));
CPoleTS = zeros(1, NUM_CTRL+1, size(x_starts, 2));
UPoleTS = zeros(length(sys_PoleTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KPoleTS = zeros(length(sys_PoleTS.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
XCartFClose = nan(4, NUM_CTRL+1, size(x_starts, 2));
UCartFClose = zeros(length(sys_CartFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCartFClose = zeros(length(sys_CartFF.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
TracePoleTS = cell(size(x_starts, 2), 1);
timePoleTS = zeros(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    dyn_PoleTS = @(x, u, k, K, xn, i) ...
               cartpole_dyn_second_cst(sys_PoleTS, x, u, k, K, xn, sys_PoleTS.full_DDP);
    Uinit = zeros(length(sys_PoleTS.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XPoleTS(:,:, jj), ...
     UPoleTS(:,1:NUM_CTRL, jj), ...
     KPoleTS(:,:,1:NUM_CTRL, jj), ...
     UCartFClose(:,:, jj), ...
     KCartFClose(:,:,:, jj), ...
     XCartFClose(:,:, jj), ...
     ~, ~, CPoleTS(:,:, jj), TracePoleTS{jj,1}] = iLQGSecondKDTree(dyn_PoleTS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UCartFF(:,:, ceil(jj/size(pole_starts, 2))), ...
                                        KCartFF(:,:,:, ceil(jj/size(pole_starts, 2))), ...
                                        XCartFF(sys_CartFF.X_DIMS_FREE,:, ceil(jj/size(pole_starts, 2))), ...
                                        Op);
    timePoleTS(jj) = toc;
    jj
end

%% Pole - T, Cart - F
disp('**** T - Pole, F - Both ****');
% Pole first
disp('T - Pole')
sys_PoleTF = sys;
sys_PoleTF.U_DIMS_FREE = [2];
sys_PoleTF.U_DIMS_FIXED = linspace(1,2,2)';
sys_PoleTF.U_DIMS_FIXED(sys_PoleTF.U_DIMS_FREE) = [];
sys_PoleTF.X_DIMS_FREE = [3; 4];
sys_PoleTF.X_DIMS_FIXED = linspace(1,4,4)';
sys_PoleTF.X_DIMS_FIXED(sys_PoleTF.X_DIMS_FREE) = [];

Op.lims = sys_PoleTF.lims(sys_PoleTF.U_DIMS_FREE, :);
XPoleTF = zeros(4, NUM_CTRL+1, size(pole_starts, 2));
CPoleTF = zeros(1, NUM_CTRL+1, size(pole_starts, 2));
UPoleTF = zeros(length(sys_PoleTF.U_DIMS_FREE), NUM_CTRL+1, size(pole_starts, 2));
KPoleTF = zeros(length(sys_PoleTF.U_DIMS_FREE), 4, NUM_CTRL+1, size(pole_starts, 2));
TracePoleTF = cell(size(pole_starts, 2), 1);
timePoleTF = zeros(size(pole_starts, 2), 1);
for jj=1:1:size(pole_starts, 2)
    dyn_PoleTF = @(x, u, i) cartpole_dyn_first_cst(sys_PoleTF, x, u, sys_PoleTF.full_DDP);
    Uinit = zeros(length(sys_PoleTF.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XPoleTF(sys_PoleTF.X_DIMS_FREE,:, jj), ...
     UPoleTF(:,1:NUM_CTRL, jj), ...
     KPoleTF(:,sys_PoleTF.X_DIMS_FREE,1:NUM_CTRL, jj), ...
     ~, ~, CPoleTF(:,:, jj), TracePoleTF{jj,1}] = iLQG(dyn_PoleTF, ...
                                    x_starts(sys_PoleTF.X_DIMS_FREE, jj), ...
                                    Uinit, Op);
    timePoleTF(jj) = toc;
    XPoleTF(sys_PoleTF.X_DIMS_FIXED,:, jj) = sys_PoleTF.l_point(sys_PoleTF.X_DIMS_FIXED)...
                                              * ones(1, size(XPoleTF, 2));
    jj
end

% Cart second
disp('F - Both')
sys_CartFS = sys;
sys_CartFS.U_DIMS_FREE = [1];
sys_CartFS.U_DIMS_FIXED = linspace(1,2,2)';
sys_CartFS.U_DIMS_FIXED(sys_CartFS.U_DIMS_FREE) = [];
sys_CartFS.X_DIMS_FREE = [1;2;3;4];
sys_CartFS.X_DIMS_FIXED = [];

Op.lims = sys_CartFS.lims(sys_CartFS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [3;4];
XCartFS = zeros(4, NUM_CTRL+1, size(x_starts, 2));
CCartFS = zeros(1, NUM_CTRL+1, size(x_starts, 2));
UCartFS = zeros(length(sys_CartFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCartFS = zeros(length(sys_CartFS.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
XPoleTClose = nan(4, NUM_CTRL+1, size(x_starts, 2));
UPoleTClose = zeros(length(sys_PoleTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KPoleTClose = zeros(length(sys_PoleTF.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
TraceCartFS = cell(size(x_starts, 2), 1);
timeCartFS = zeros(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    dyn_CartFS = @(x, u, k, K, xn, i) ...
               cartpole_dyn_second_cst(sys_CartFS, x, u, k, K, xn, sys_CartFS.full_DDP);
    Uinit = zeros(length(sys_CartFS.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XCartFS(:,:, jj), ...
     UCartFS(:,1:NUM_CTRL, jj), ...
     KCartFS(:,:,1:NUM_CTRL, jj), ...
     UPoleTClose(:,:, jj), ...
     KPoleTClose(:,:,:, jj), ...
     XPoleTClose(:,:, jj), ...
     ~, ~, CCartFS(:,:, jj), TraceCartFS{jj,1}] = iLQGSecondKDTree(dyn_CartFS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UPoleTF(:,:, mod(jj-1, size(pole_starts,2))+1), ...
                                        KPoleTF(:,:,:, mod(jj-1, size(pole_starts,2))+1), ...
                                        XPoleTF(sys_PoleTF.X_DIMS_FREE,:, mod(jj-1, size(pole_starts,2))+1), ...
                                        Op);
    timeCartFS(jj) = toc;
    jj
end

%% Cart - T, Pole - F
disp('**** T - Cart, F - Both ****');
% Cart first
disp('T - Cart')
sys_CartTF = sys;
sys_CartTF.U_DIMS_FREE = [2];
sys_CartTF.U_DIMS_FIXED = linspace(1,2,2)';
sys_CartTF.U_DIMS_FIXED(sys_CartTF.U_DIMS_FREE) = [];
sys_CartTF.X_DIMS_FREE = [1;2];
sys_CartTF.X_DIMS_FIXED = linspace(1,4,4)';
sys_CartTF.X_DIMS_FIXED(sys_CartTF.X_DIMS_FREE) = [];

Op.lims = sys_CartTF.lims(sys_CartTF.U_DIMS_FREE, :);
XCartTF = zeros(4, NUM_CTRL+1, size(cart_starts, 2));
CCartTF = zeros(1, NUM_CTRL+1, size(cart_starts, 2));
UCartTF = zeros(length(sys_CartTF.U_DIMS_FREE), NUM_CTRL+1, size(cart_starts, 2));
KCartTF = zeros(length(sys_CartTF.U_DIMS_FREE), 4, NUM_CTRL+1, size(cart_starts, 2));
TraceCartTF = cell(size(cart_starts, 2), 1);
timeCartTF = zeros(size(cart_starts, 2), 1);
for jj=1:1:size(cart_starts, 2)
    dyn_CartTF = @(x, u, i) cartpole_dyn_first_cst(sys_CartTF, x, u, sys_CartTF.full_DDP);
    Uinit = zeros(length(sys_CartTF.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XCartTF(sys_CartTF.X_DIMS_FREE,:, jj), ...
     UCartTF(:,1:NUM_CTRL, jj), ...
     KCartTF(:,sys_CartTF.X_DIMS_FREE,1:NUM_CTRL, jj), ...
     ~, ~, CCartTF(:,:, jj), TraceCartTF{jj,1}] = iLQG(dyn_CartTF, ...
                                    x_starts(sys_CartTF.X_DIMS_FREE, jj), ...
                                    Uinit, Op);
    timeCartTF(jj) = toc;
    XCartTF(sys_CartTF.X_DIMS_FIXED,:, jj) = sys_CartTF.l_point(sys_CartTF.X_DIMS_FIXED)...
                                              * ones(1, size(XCartTF, 2));
    jj
end

% Pole second
disp('F - Both')
sys_PoleFS = sys;
sys_PoleFS.U_DIMS_FREE = [1];
sys_PoleFS.U_DIMS_FIXED = linspace(1,2,2)';
sys_PoleFS.U_DIMS_FIXED(sys_PoleFS.U_DIMS_FREE) = [];
sys_PoleFS.X_DIMS_FREE = [1;2;3;4];
sys_PoleFS.X_DIMS_FIXED = [];

Op.lims = sys_PoleFS.lims(sys_PoleFS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [1;2];
XPoleFS = zeros(4, NUM_CTRL+1, size(x_starts, 2));
CPoleFS = zeros(1, NUM_CTRL+1, size(x_starts, 2));
UPoleFS = zeros(length(sys_PoleFS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KPoleFS = zeros(length(sys_PoleFS.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
XCartTClose = nan(4, NUM_CTRL+1, size(x_starts, 2));
UCartTClose = zeros(length(sys_CartTF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCartTClose = zeros(length(sys_CartTF.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
TracePoleFS = cell(size(x_starts, 2), 1);
timePoleFS = zeros(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    dyn_PoleFS = @(x, u, k, K, xn, i) ...
               cartpole_dyn_second_cst(sys_PoleFS, x, u, k, K, xn, sys_PoleFS.full_DDP);
    Uinit = zeros(length(sys_PoleFS.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XPoleFS(:,:, jj), ...
     UPoleFS(:,1:NUM_CTRL, jj), ...
     KPoleFS(:,:,1:NUM_CTRL, jj), ...
     UCartTClose(:,:, jj), ...
     KCartTClose(:,:,:, jj), ...
     XCartTClose(:,:, jj), ...
     ~, ~, CPoleFS(:,:, jj), TracePoleFS{jj,1}] = iLQGSecondKDTree(dyn_PoleFS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UCartTF(:,:, ceil(jj/size(pole_starts, 2))), ...
                                        KCartTF(:,:,:, ceil(jj/size(pole_starts, 2))), ...
                                        XCartTF(sys_CartTF.X_DIMS_FREE,:, ceil(jj/size(pole_starts, 2))), ...
                                        Op);
    timePoleFS(jj) = toc;
    jj
end

%% Pole - F, Cart - T
disp('**** F - Pole, T - Both ****');
% Pole first
disp('F - Pole')
sys_PoleFF = sys;
sys_PoleFF.U_DIMS_FREE = [1];
sys_PoleFF.U_DIMS_FIXED = linspace(1,2,2)';
sys_PoleFF.U_DIMS_FIXED(sys_PoleFF.U_DIMS_FREE) = [];
sys_PoleFF.X_DIMS_FREE = [3; 4];
sys_PoleFF.X_DIMS_FIXED = linspace(1,4,4)';
sys_PoleFF.X_DIMS_FIXED(sys_PoleFF.X_DIMS_FREE) = [];

Op.lims = sys_PoleFF.lims(sys_PoleFF.U_DIMS_FREE, :);
XPoleFF = zeros(4, NUM_CTRL+1, size(pole_starts, 2));
CPoleFF = zeros(1, NUM_CTRL+1, size(pole_starts, 2));
UPoleFF = zeros(length(sys_PoleFF.U_DIMS_FREE), NUM_CTRL+1, size(pole_starts, 2));
KPoleFF = zeros(length(sys_PoleFF.U_DIMS_FREE), 4, NUM_CTRL+1, size(pole_starts, 2));
TracePoleFF = cell(size(pole_starts, 2), 1);
timePoleFF = zeros(size(pole_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    dyn_PoleFF = @(x, u, i) cartpole_dyn_first_cst(sys_PoleFF, x, u, sys_PoleFF.full_DDP);
    Uinit = zeros(length(sys_PoleFF.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XPoleFF(sys_PoleFF.X_DIMS_FREE,:, jj), ...
     UPoleFF(:,1:NUM_CTRL, jj), ...
     KPoleFF(:,sys_PoleFF.X_DIMS_FREE,1:NUM_CTRL, jj), ...
     ~, ~, CPoleFF(:,:, jj), TracePoleFF{jj,1}] = iLQG(dyn_PoleFF, ...
                                                    x_starts(sys_PoleFF.X_DIMS_FREE, jj), ...
                                                    Uinit, Op);
    timePoleFF(jj) = toc;
    XPoleFF(sys_PoleFF.X_DIMS_FIXED,:, jj) = sys_PoleFF.l_point(sys_PoleFF.X_DIMS_FIXED)...
                                              * ones(1, size(XPoleFF, 2));
    jj
end

% Cart second
disp('T - Both')
sys_CartTS = sys;
sys_CartTS.U_DIMS_FREE = [2];
sys_CartTS.U_DIMS_FIXED = linspace(1,2,2)';
sys_CartTS.U_DIMS_FIXED(sys_CartTS.U_DIMS_FREE) = [];
sys_CartTS.X_DIMS_FREE = [1;2;3;4];
sys_CartTS.X_DIMS_FIXED = [];

Op.lims = sys_CartTS.lims(sys_CartTS.U_DIMS_FREE, :);
Op.X_DIMS_FIRST = [3;4];
XCartTS = zeros(4, NUM_CTRL+1, size(x_starts, 2));
CCartTS = zeros(1, NUM_CTRL+1, size(x_starts, 2));
UCartTS = zeros(length(sys_CartTS.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KCartTS = zeros(length(sys_CartTS.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
XPoleFClose = nan(4, NUM_CTRL+1, size(x_starts, 2));
UPoleFClose = zeros(length(sys_PoleFF.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
KPoleFClose = zeros(length(sys_PoleFF.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
TraceCartTS = cell(size(x_starts, 2), 1);
timeCartTS = zeros(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    dyn_CartTS = @(x, u, k, K, xn, i) ...
               cartpole_dyn_second_cst(sys_CartTS, x, u, k, K, xn, sys_CartTS.full_DDP);
    Uinit = zeros(length(sys_CartTS.U_DIMS_FREE), NUM_CTRL);
    tic;
    [XCartTS(:,:, jj), ...
     UCartTS(:,1:NUM_CTRL, jj), ...
     KCartTS(:,:,1:NUM_CTRL, jj), ...
     UPoleFClose(:,:, jj), ...
     KPoleFClose(:,:,:, jj), ...
     XPoleFClose(:,:, jj), ...
     ~, ~, CCartTS(:,:, jj), TraceCartTS{jj,1}] = iLQGSecondKDTree(dyn_CartTS, ...
                                        x_starts(:, jj), ...
                                        Uinit, ...
                                        UPoleFF(:,:, ceil(jj/size(pole_starts, 2))), ...
                                        KPoleFF(:,:,:, ceil(jj/size(pole_starts, 2))), ...
                                        XPoleFF(sys_PoleFF.X_DIMS_FREE,:, ceil(jj/size(pole_starts, 2))), ...
                                        Op);
    timeCartTS(jj) = toc;
    jj
end

%% Decoupled
% Cart - F, Pole - T
disp('**** F - Cart, T - Pole ****');
sys_CartFPoleTDec = sys;
sys_CartFPoleTDec.U_DIMS_FREE = [1;2];
sys_CartFPoleTDec.U_DIMS_FIXED = [];
sys_CartFPoleTDec.X_DIMS_FREE = [1;2;3;4];
sys_CartFPoleTDec.X_DIMS_FIXED = [];
XCartFPoleTDec = zeros(4, NUM_CTRL+1, size(x_starts, 2));
UCartFPoleTDec = zeros(2, NUM_CTRL, size(x_starts, 2));
timeCartFPoleT = zeros(size(x_starts, 2), 1);
U_DIMS = cell(2,1);
U_DIMS{1,1} = [1];
U_DIMS{2,1} = [2];

for jj=1:1:size(x_starts, 2)
    tic;
    dyn_CartFPoleTDec = @(x, u) ...
               cartpole_dyn_first_cst(sys_CartFPoleTDec, x, u, sys_CartFPoleTDec.full_DDP);
    [XCartFPoleTDec(:,:, jj), UCartFPoleTDec(:,:, jj)] = ForwardPassDec([UCartFF(:,:, ceil(jj/size(pole_starts, 2))); UPoleTF(:,:, mod(jj-1, size(pole_starts,2))+1)], ...
                                                                        [KCartFF(:,:,:, ceil(jj/size(pole_starts, 2))); KPoleTF(:,:,:, mod(jj-1, size(pole_starts,2))+1)], ...
                                                                        [XCartFF(:,:, ceil(jj/size(pole_starts, 2))); XPoleTF(:,:, mod(jj-1, size(pole_starts,2))+1)], ...
                                                                        U_DIMS, ... 
                                                                        x_starts(:, jj), dyn_CartFPoleTDec);
    timeCartFPoleT(jj) = toc;
    jj
end

% Cart - T, Pole - F
disp('**** T - Cart, F - Pole ****');
sys_CartTPoleFDec = sys;
sys_CartTPoleFDec.U_DIMS_FREE = [1;2];
sys_CartTPoleFDec.U_DIMS_FIXED = [];
sys_CartTPoleFDec.X_DIMS_FREE = [1;2;3;4];
sys_CartTPoleFDec.X_DIMS_FIXED = [];
XCartTPoleFDec = zeros(4, NUM_CTRL+1, size(x_starts, 2));
UCartTPoleFDec = zeros(2, NUM_CTRL, size(x_starts, 2));
timeCartTPoleF = zeros(size(x_starts, 2), 1);
for jj=1:1:size(x_starts, 2)
    tic;
    dyn_CartTPoleFDec = @(x, u) ...
               cartpole_dyn_first_cst(sys_CartTPoleFDec, x, u, sys_CartTPoleFDec.full_DDP);
    [XCartTPoleFDec(:,:, jj), UCartTPoleFDec(:,:, jj)] = ForwardPassDec([UPoleFF(:,:, mod(jj-1, size(pole_starts,2))+1); UCartTF(:,:, ceil(jj/size(pole_starts, 2)))], ...
                                                                        [KPoleFF(:,:,:, mod(jj-1, size(pole_starts,2))+1); KCartTF(:,:,:, ceil(jj/size(pole_starts, 2)))], ...
                                                                        [XPoleFF(:,:, mod(jj-1, size(pole_starts,2))+1); XCartTF(:,:, ceil(jj/size(pole_starts, 2)))], ...
                                                                        U_DIMS, ...   
                                                                        x_starts(:, jj), dyn_CartTPoleFDec);
    timeCartTPoleF(jj) = toc;
    jj
end

status = 'zeroinit_finallpoint_nowraparound_difffurthercloserstarts_finestalpha_coarsergamma';

save(strcat('data/iLQGCartPoleDecomposed_diffR_difffurthercloserstarts_alpha_',num2str(Op.Alpha),'_gamma_',num2str(sys.gamma_),'mc=',num2str(sys.mc),',mp=',num2str(sys.mp),'.mat'));