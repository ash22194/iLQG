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
Op.Alpha = [1];

% Define starts
%cart_starts = [-1, -0.5, 0, 0.5, 1;
%                0,  0,   0, 0,   0];
%cart_starts = [-0.75, -0.5, 0, 0.5, 0.75;
%                0,     0,   0, 0,   0];
cart_starts = [-0.4, -0.2, 0, 0.2, 0.4;
                0,     0,   0, 0,   0];

%pole_starts = [7*pi/4, 5*pi/4, 3*pi/4, pi/4, 0;
%               0,      0,      0,      0,    0];
%pole_starts = [pi/2, 2*pi/3, 5*pi/6, 7*pi/6, 4*pi/3, 3*pi/2;
%               0,      0,      0,      0,    0,      0];
pole_starts = [2*pi/3, 3*pi/4, 5*pi/6, 7*pi/6, 5*pi/4, 4*pi/3;
               0,      0,      0,      0,    0,      0];

x_starts = nan(4, size(cart_starts,2)*size(pole_starts,2));
for ii=1:1:size(cart_starts, 2)
    for jj=1:1:size(pole_starts, 2)
        x_starts(:, (ii-1)*size(pole_starts, 2) + jj) = [cart_starts(:, ii); pole_starts(:, jj)];
    end
end

x_starts = [0; 0; 0; 0];

%% Cascaded

X_DIMS = [1;2;3;4];
U_DIMS = [1;2];

sysFF = cell(2^length(X_DIMS) - 2, 1);
sysTS = cell(2^length(X_DIMS) - 2, 1);

% F First
disp("Cascaded F First");
numCompleted = 0;
for kk=1:3
    C = nchoosek(X_DIMS, kk);
    for ii=1:1:size(C,1)
        sysFF_ = sys;
        sysFF_.X_DIMS_FREE = C(ii,:);
        sysFF_.X_DIMS_FIXED = X_DIMS;
        sysFF_.X_DIMS_FIXED(sysFF_.X_DIMS_FREE) = [];
        sysFF_.U_DIMS_FREE = [1];
        sysFF_.U_DIMS_FIXED = [2];
        Op.lims = sysFF_.lims(sysFF_.U_DIMS_FREE, :);

        sysFF_.Xn = zeros(4, NUM_CTRL+1, size(x_starts, 2));
        sysFF_.Cn = zeros(1, NUM_CTRL+1, size(x_starts, 2));
        sysFF_.Un = zeros(length(sysFF_.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
        sysFF_.Kn = zeros(length(sysFF_.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
        sysFF_.timen = zeros(size(x_starts, 2), 1);

        dynFF = @(x, u, i) cartpole_dyn_first_cst(sysFF_, x, u, sysFF_.full_DDP);
        for jj=1:1:size(x_starts, 2)

            tic;
            [sysFF_.Xn(sysFF_.X_DIMS_FREE, :, jj), ...
             sysFF_.Un(:, 1:NUM_CTRL, jj), ...
             sysFF_.Kn(:, sysFF_.X_DIMS_FREE, 1:NUM_CTRL, jj), ...
             ~, ~, sysFF_.Cn(:, :, jj)] = iLQG(dynFF, ...
                                            x_starts(sysFF_.X_DIMS_FREE, jj), ...
                                            zeros(length(sysFF_.U_DIMS_FREE), NUM_CTRL), Op);
            sysFF_.Xn(sysFF_.X_DIMS_FIXED,:, jj) = sysFF_.l_point(sysFF_.X_DIMS_FIXED) ...
                                                      * ones(1, NUM_CTRL+1);
            sysFF_.timen(jj) = toc;
            disp(strcat("FF Num dims :", num2str(kk), ", Trajectory :", num2str(jj)));
        end
        sysFF{numCompleted + ii, 1} = sysFF_;

        sysTS_ = sys;
        sysTS_.X_DIMS_FREE = X_DIMS;
        sysTS_.X_DIMS_FIXED = [];
        sysTS_.U_DIMS_FREE = [2];
        sysTS_.U_DIMS_FIXED = [1];
        Op.lims = sysTS_.lims(sysTS_.U_DIMS_FREE, :);
        Op.X_DIMS_FIRST = sysFF_.X_DIMS_FREE;

        sysTS_.Xn = zeros(4, NUM_CTRL+1, size(x_starts, 2));
        sysTS_.Cn = zeros(1, NUM_CTRL+1, size(x_starts, 2));
        sysTS_.Un = zeros(length(sysTS_.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
        sysTS_.Kn = zeros(length(sysTS_.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
        sysTS_.XFC = nan(4, NUM_CTRL+1, size(x_starts, 2));
        sysTS_.UFC = zeros(length(sysFF_.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
        sysTS_.KFC = zeros(length(sysFF_.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
        sysTS_.timen = zeros(size(x_starts, 2), 1);

        dynTS = @(x, u, k, K, xn, i) ...
                   cartpole_dyn_second_cst(sysTS_, x, u, k, K, xn, sysTS_.full_DDP);
        for jj=1:1:size(x_starts, 2)

            tic;
            [sysTS_.Xn(:,:, jj), ...
             sysTS_.Un(:,1:NUM_CTRL, jj), ...
             sysTS_.Kn(:,:,1:NUM_CTRL, jj), ...
             sysTS_.XFC(:,:, jj), ...
             sysTS_.UFC(:,:,:, jj), ...
             sysTS_.KFC(:,:, jj), ...
             ~, ~, sysTS_.Cn(:,:, jj)] = iLQGSecondKDTree(dynTS, ...
                                                x_starts(:, jj), ...
                                                zeros(length(sysTS_.U_DIMS_FREE), NUM_CTRL), ...
                                                sysFF_.Un(:,:, jj), ...
                                                sysFF_.Kn(:,:,:, jj), ...
                                                sysFF_.Xn(sysFF_.X_DIMS_FREE,:, jj), ...
                                                Op);
            sysTS_.timen(jj) = toc;
            disp(strcat("TS Num dims :", num2str(4 - kk), ", Trajectory :", num2str(jj)));
        end
        sysTS{numCompleted + ii, 1} = sysTS_;
    end
    numCompleted = numCompleted + size(C, 1);
end

sysTF = cell(2^length(X_DIMS) - 2, 1);
sysFS = cell(2^length(X_DIMS) - 2, 1);

% T First
disp("Cascaded T First");
numCompleted = 0;
for kk=1:3
    C = nchoosek(X_DIMS, kk);
    for ii=1:1:size(C,1)
        sysTF_ = sys;
        sysTF_.X_DIMS_FREE = C(ii,:);
        sysTF_.X_DIMS_FIXED = X_DIMS;
        sysTF_.X_DIMS_FIXED(sysTF_.X_DIMS_FREE) = [];
        sysTF_.U_DIMS_FREE = [2];
        sysTF_.U_DIMS_FIXED = [1];
        Op.lims = sysTF_.lims(sysTF_.U_DIMS_FREE, :);

        sysTF_.Xn = zeros(4, NUM_CTRL+1, size(x_starts, 2));
        sysTF_.Cn = zeros(1, NUM_CTRL+1, size(x_starts, 2));
        sysTF_.Un = zeros(length(sysTF_.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
        sysTF_.Kn = zeros(length(sysTF_.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
        sysTF_.timen = zeros(size(x_starts, 2), 1);

        dynTF = @(x, u, i) cartpole_dyn_first_cst(sysTF_, x, u, sysTF_.full_DDP);
        for jj=1:1:size(x_starts, 2)

            tic;
            [sysTF_.Xn(sysTF_.X_DIMS_FREE, :, jj), ...
             sysTF_.Un(:, 1:NUM_CTRL, jj), ...
             sysTF_.Kn(:, sysTF_.X_DIMS_FREE, 1:NUM_CTRL, jj), ...
             ~, ~, sysTF_.Cn(:, :, jj)] = iLQG(dynTF, ...
                                            x_starts(sysTF_.X_DIMS_FREE, jj), ...
                                            zeros(length(sysTF_.U_DIMS_FREE), NUM_CTRL), Op);
            sysTF_.Xn(sysTF_.X_DIMS_FIXED,:, jj) = sysTF_.l_point(sysTF_.X_DIMS_FIXED) ...
                                                      * ones(1, NUM_CTRL+1);
            sysTF_.timen(jj) = toc;
            disp(strcat("TF Num dims :", num2str(kk), ", Trajectory :", num2str(jj)));
        end
        sysTF{numCompleted + ii, 1} = sysTF_;

        sysFS_ = sys;
        sysFS_.X_DIMS_FREE = X_DIMS;
        sysFS_.X_DIMS_FIXED = [];
        sysFS_.U_DIMS_FREE = [1];
        sysFS_.U_DIMS_FIXED = [2];
        Op.lims = sysFS_.lims(sysFS_.U_DIMS_FREE, :);
        Op.X_DIMS_FIRST = sysTF_.X_DIMS_FREE;

        sysFS_.Xn = zeros(4, NUM_CTRL+1, size(x_starts, 2));
        sysFS_.Cn = zeros(1, NUM_CTRL+1, size(x_starts, 2));
        sysFS_.Un = zeros(length(sysFS_.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
        sysFS_.Kn = zeros(length(sysFS_.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
        sysFS_.XFC = nan(4, NUM_CTRL+1, size(x_starts, 2));
        sysFS_.UFC = zeros(length(sysTF_.U_DIMS_FREE), NUM_CTRL+1, size(x_starts, 2));
        sysFS_.KFC = zeros(length(sysTF_.U_DIMS_FREE), 4, NUM_CTRL+1, size(x_starts, 2));
        sysFS_.timen = zeros(size(x_starts, 2), 1);

        dynFS = @(x, u, k, K, xn, i) ...
                   cartpole_dyn_second_cst(sysFS_, x, u, k, K, xn, sysFS_.full_DDP);
        for jj=1:1:size(x_starts, 2)

            tic;
            [sysFS_.Xn(:,:, jj), ...
             sysFS_.Un(:,1:NUM_CTRL, jj), ...
             sysFS_.Kn(:,:,1:NUM_CTRL, jj), ...
             sysFS_.XFC(:,:, jj), ...
             sysFS_.UFC(:,:,:, jj), ...
             sysFS_.KFC(:,:, jj), ...
             ~, ~, sysFS_.Cn(:,:, jj)] = iLQGSecondKDTree(dynFS, ...
                                                x_starts(:, jj), ...
                                                zeros(length(sysFS_.U_DIMS_FREE), NUM_CTRL), ...
                                                sysTF_.Un(:,:, jj), ...
                                                sysTF_.Kn(:,:,:, jj), ...
                                                sysTF_.Xn(sysTF_.X_DIMS_FREE,:, jj), ...
                                                Op);
            sysFS_.timen(jj) = toc;
            disp(strcat("FS Num dims :", num2str(4 - kk), ", Trajectory :", num2str(jj)));
        end
        sysFS{numCompleted + ii, 1} = sysFS_;
    end
    numCompleted = numCompleted + size(C, 1);
end
