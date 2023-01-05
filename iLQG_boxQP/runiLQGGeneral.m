clear;
close all;
clc;

%% 

system_name = 'cartpole';
restoredefaultpath();
addpath(strcat('new_systems/', system_name));
addpath('general');

decompositions = load(strcat('../GA/data/', system_name, 'AllPossibleDecompositions.mat'));
sys = decompositions.sys;
sys.name = strcat(system_name, '_pareto');
sys.cxmod = zeros(sys.X_DIMS,1);
sys.full_DDP = false;
sys.u0init = false;

Op.lims  = sys.lims;
Op.maxIter = 500;
Op.gamma_ = sys.gamma_;
Op.print = 0;
Op.save_dir = 'data/multitraj';
Op.reuse_policy = true;

if (strcmp(system_name, 'cartpole'))
    cart_starts = [-0.5, -0.5, 0.5, 0.5;
                   -1,  1, -1, 1];
    cart_starts(2,:) = 0.5 * cart_starts(2,:);
    pole_starts = [2*pi/3, 2*pi/3, 4*pi/3, 4*pi/3;
                   -1,    1,     -1,   1];
    pole_starts(2,:) = 0.5 * pole_starts(2,:);
    starts = zeros(sys.X_DIMS, size(cart_starts, 2) ...
                               *size(pole_starts, 2));
    count = 0;
    for carts=1:1:size(cart_starts, 2)
        for poles=1:1:size(pole_starts, 2)
            count = count + 1;
            starts(1:2, count) = cart_starts(:, carts);
            starts(3:4, count) = pole_starts(:, poles);
        end
    end
    sys.T = 4;
    
elseif (strcmp(system_name, 'manipulator2dof') || strcmp(system_name, 'manipulator3dof'))
    
    if (strcmp(system_name, 'manipulator2dof'))
        n = 2;
    elseif (strcmp(system_name, 'manipulator3dof'))
        n = 3;
    end
    
    theta_starts(:,:,1) = [2*pi/3, 2*pi/3, 4*pi/3, 4*pi/3;
                           -1,     1,      -1,     1];
    theta_starts(:,:,2:n) = repmat([-pi/3, -pi/3, pi/3, pi/3;
                                    -1,    1,     -1,   1], [1,1,n-1]);
    theta_starts(2,:,:) = 0.5*theta_starts(2,:,:);

    starts = zeros(sys.X_DIMS, size(theta_starts, 2)^n);
    for count = 1:1:size(starts,2)
        count_ = count;
        for jj=1:1:n
            index = 1 + mod(count_-1, size(theta_starts, 2));
            starts([jj; n+jj], count) = theta_starts(:, index, jj);
            count_ = (count_ - index)/size(theta_starts, 2);
        end
    end
    sys.T = 4;

elseif (strcmp(system_name, 'biped2d'))
    com_pos = [0.92, 0.92, 1.0, 1.0;
               0.3,  0.2, 0.3, 0.2]; % adjusted state bounds
    com_pos(2,:) = pi/2 + com_pos(2,:);
    com_vel = [ 0.1, -0.1, 0.1, -0.1;
               -0.3, -0.3, 0.3, 0.3];
    theta_starts = [-0.2,  -0.2, 0.2,  0.2;
                    -0.2,   0.2, -0.2, 0.2];

    starts = zeros(sys.X_DIMS, size(com_pos, 2) ...
                               *size(com_vel, 2) ...
                               *size(theta_starts, 2));
    count = 0;
    for cp=1:1:size(com_pos, 2)
        for cv=1:1:size(com_vel, 2)
            for ts=1:1:size(theta_starts, 2)
                count = count + 1;
                starts(1:2, count) = com_pos(:, cp);
                starts(3:4, count) = com_vel(:, cv);
                starts(5:6, count) = theta_starts(:, ts);
            end
        end
    end
    sys.T = 4;
    
end

if (size(decompositions.policy_decompositions, 1) > 44)
    policy_decompositions = decompositions.policy_decompositions_pareto_front;
else
    policy_decompositions = decompositions.policy_decompositions;
end

Xd = zeros(sys.X_DIMS, round(sys.T / sys.dt)+1, size(starts, 2), size(policy_decompositions, 1));
Ud = zeros(sys.U_DIMS, round(sys.T / sys.dt), size(starts, 2), size(policy_decompositions, 1));
cd = zeros(1, size(starts, 2), size(policy_decompositions, 1));

for dd=1:1:size(policy_decompositions, 1)
    p = policy_decompositions{dd,1};
    s = policy_decompositions{dd,2};
    fprintf('Decomposition : %d / %d\n', dd, size(policy_decompositions, 1));
    sys.decomposition_id = dd;
    [Xd(:,:,:,dd), Ud(:,:,:,dd), cd(:,:,dd)] = ilqg_decomposition_multitraj(sys, Op, p, s, starts);
end

p_joint = [zeros(sys.U_DIMS,1), ones(sys.U_DIMS,1)];
s_joint = ones(sys.U_DIMS, sys.X_DIMS);
sys.decomposition_id = 0;
[Xjoint, Ujoint, cjoint] = ilqg_decomposition_multitraj(sys, Op, p_joint, s_joint, starts);

err_ddp = mean(reshape(cd, size(starts, 2), size(policy_decompositions, 1)) - cjoint', 1);
