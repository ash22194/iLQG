clear;
clc;
close all;

%%
addpath('cartpole');
iLQG_dir = 'data/';
iLQG_AllPossible_filename = "iLQGCartPoleAllPossibleDecomposed_diffR_furthercloserstarts_alpha_0.1_gamma_0.997mc=5,mp=1.mat";
iLQG_Joint_filename = "iLQGCartPoleAllPossibleDecomposed_diffR_furthercloserstarts_alpha_0.1_gamma_0.997mc=5,mp=1.mat";

load(strcat(iLQG_dir, iLQG_AllPossible_filename), 'sysFF', 'sysTS', 'sysTF', 'sysFS', 'x_starts', 'NUM_CTRL');
load(strcat(iLQG_dir, iLQG_Joint_filename), 'sysJoint', 'sys');
CJoint = sysJoint.Cn;

%% Cascaded
% include_trajectories = [1, 2, 5, 6, 7, 8, 11, 12, 19, 20, 23, 24, 25, 26, 29, 30];
include_trajectories = linspace(1,16,16);
JJoint = reshape(sum(CJoint(1,:,:), 2), length(include_trajectories), 1);

% F First
CascFF = cell(size(sysTS, 1), 1);
JErrFF = nan(size(sysTS, 1), 1);

for ii=1:1:size(sysTS, 1)
    sysFF_ = sysFF{ii,1};
    sysTS_ = sysTS{ii,1};
    JTS = reshape(sum(sysTS_.Cn(1,:,:), 2), size(x_starts, 2), 1);
    CascFF{ii, 1} = sysFF_.X_DIMS_FREE;
    JErrFF(ii, 1) = sum(abs(JTS(include_trajectories, 1) - JJoint));
    if (isempty(setdiff(sysFF_.X_DIMS_FREE, [3,4])) && isempty(setdiff([3,4], sysFF_.X_DIMS_FREE)))
       disp('wait'); 
    end
    ii
end

% T First
CascTF = cell(size(sysFS, 1), 1);
JErrTF = nan(size(sysFS, 1), 1);
for ii=1:1:size(sysFS, 1)
    sysTF_ = sysTF{ii,1};
    sysFS_ = sysFS{ii,1};
    JFS = reshape(sum(sysFS_.Cn(1,:,:), 2), size(x_starts, 2), 1);
    CascTF{ii, 1} = sysTF_.X_DIMS_FREE;
    JErrTF(ii, 1) = sum(abs(JFS(include_trajectories, 1) - JJoint));
    ii
end

%% Decoupled

Dec = cell(size(sysTS, 1)-1, 1);
JErrDec = nan(size(sysTS, 1)-1, 1);
sys.X_DIMS_FREE = [1;2;3;4];
sys.X_DIMS_FIXED = [];
sys.U_DIMS_FREE = [1;2];
sys.U_DIMS_FIXED = [];

U_DIMS = cell(2,1);
U_DIMS{1,1} = [1];
U_DIMS{2,1} = [2];
for ii=1:1:(size(sysTS, 1)-1)
    sysFF_ = sysFF{ii,1};
    F_DIMS_FREE = sysFF_.X_DIMS_FREE;
    Dec{ii, 1} = F_DIMS_FREE;
    T_DIMS_FREE = linspace(1,4,4);
    T_DIMS_FREE(F_DIMS_FREE) = [];
    
    for jj=1:1:size(sysTF, 1)
        sysTF_ = sysTF{jj,1};
        if (isempty(setdiff(T_DIMS_FREE, sysTF_.X_DIMS_FREE)) ...
                && isempty(setdiff(sysTF_.X_DIMS_FREE, T_DIMS_FREE)))
            break;
        end
    end
    
    X = zeros(4, NUM_CTRL+1, length(include_trajectories));
    U = zeros(2, NUM_CTRL, length(include_trajectories));
    JDec = zeros(length(include_trajectories), 1);
    count = 1;
    for jj=include_trajectories
        
        dyn_Dec = @(x, u) ...
                   cartpole_dyn_first_cst(sys, x, u, sys.full_DDP);
        [X(:,:, count), U(:,:, count)] = ForwardPassDec([sysFF_.Un(:,:, jj); sysTF_.Un(:,:, jj)], ...
                                                  [sysFF_.Kn(:,:,:, jj); sysTF_.Kn(:,:,:, jj)], ...
                                                  [sysFF_.Xn(:,:, jj); sysTF_.Xn(:,:, jj)], ...
                                                   U_DIMS, ... 
                                                   x_starts(:, jj), dyn_Dec);
        JDec(count) = J_CartPole(sys, X(:,:, count), U(:,:, count));
        count
        count = count + 1;
    end
    ii
    JErrDec(ii, 1) = sum(abs(JDec - JJoint));
end

Casc = cell(size(CascFF,1) + size(CascTF,1) + size(Dec,1), 2);
Casc(1:size(CascFF,1), 1) = CascFF;
Casc(1:size(CascFF,1), 2) = {1};
Casc((size(CascFF,1) + 1):(size(CascFF,1) + size(CascTF,1)), 1) = CascTF;
Casc((size(CascFF,1) + 1):(size(CascFF,1) + size(CascTF,1)), 2) = {2};
Casc((size(CascFF,1) + size(CascTF,1) + 1):(size(CascFF,1) + size(CascTF,1) + size(Dec,1)), 1) = Dec;
Casc((size(CascFF,1) + size(CascTF,1) + 1):(size(CascFF,1) + size(CascTF,1) + size(Dec,1)), 2) = {3};
JErr = [JErrFF; JErrTF; JErrDec];
[JErrsorted, Errorder] = sort(JErr);

for nn=1:1:size(JErr, 1)
    decomp = Casc{Errorder(nn),1};
    decomptype = Casc{Errorder(nn),2};
    disp(strcat(num2str(decomp'), ', (', num2str(decomptype),'), Err : ', num2str(JErrsorted(nn))));
end
