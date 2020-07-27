clear;
clc;
close all;
addpath('cartpole');

%% Linear Dynamics

sys.mc = 5; 
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81;
sys.dt = 0.001;

% Continuous time dynamics
x = sym('x', [4, 1]);
u = sym('u', [2, 1]);
zero_dyn = [x(2);
            (sys.mp*sys.l*(x(4)^2)*sin(x(3)) + 0.5*sys.mp*sys.g*sin(2*x(3)))/(sys.mc + sys.mp*sin(x(3))^2);
            x(4);
            (-0.5*sys.mp*sys.l*x(4)^2*sin(2*x(3)) - (sys.mc+sys.mp)*sys.g*sin(x(3)))/(sys.l*(sys.mc + sys.mp*sin(x(3))^2))];
act_dyn = [0;
           (u(1) - u(2)*cos(x(3))/sys.l)/(sys.mc + sys.mp*sin(x(3))^2);
           0;
           (u(2)/sys.l*(sys.mc/sys.mp+1) - u(1)*cos(x(3)))/(sys.l*(sys.mc + sys.mp*sin(x(3))^2))];
zero_dyn_x = jacobian(zero_dyn,x);
act_dyn_u = jacobian(act_dyn,u);
sys.l_point = [0; 0; pi; 0];
sys.goal = sys.l_point;
A = eval(subs(zero_dyn_x, x, sys.l_point)); 
B = eval(subs(act_dyn_u, x, sys.l_point));

Q = diag([25, 0.02, 25, 0.02]);
R = 0.001*eye(2);
sys.gamma_ = 0.997;
sys.lambda_ = (1 - sys.gamma_)/sys.dt;

sys.numPoints = [31, 31, 31, 31];
sys.limits = [-1.5, 1.5;
              -3, 3;
              0, 2*pi;
              -3, 3];
[grid_x, grid_x_dot, grid_xP, grid_xP_dot] = ndgrid(linspace(sys.limits(1,1), sys.limits(1,2), sys.numPoints(1)), ...
                                                    linspace(sys.limits(2,1), sys.limits(2,2), sys.numPoints(2)), ...
                                                    linspace(sys.limits(3,1), sys.limits(3,2), sys.numPoints(3)), ...
                                                    linspace(sys.limits(4,1), sys.limits(4,2), sys.numPoints(4)));
sys.u_limits = [-9, 9;
                -9, 9];

max_iter = 3000;
gtol = 0.00001;

% Joint
sys_joint = sys;
sys_joint.A = A;
sys_joint.B = B;
sys_joint.Q = Q;
sys_joint.R = R;

[sys_joint.K, sys_joint.S, sys_joint.e] = lqr(A - sys.lambda_*eye(size(A,1))/2, B, Q, R, ...
                                              zeros(size(A,1), size(B,2)));
policyF = -(sys_joint.K(1,1)*(grid_x - sys.goal(1)) + sys_joint.K(1,2)*(grid_x_dot - sys.goal(2)) ...
                    + sys_joint.K(1,3)*(grid_xP - sys.goal(3)) + sys_joint.K(1,4)*(grid_xP_dot - sys.goal(4)));
policyF_bounded = min(sys.u_limits(1,2), max(sys.u_limits(1,1), policyF));
policyT = -(sys_joint.K(2,1)*(grid_x - sys.goal(1)) + sys_joint.K(2,2)*(grid_x_dot - sys.goal(2)) ...
                    + sys_joint.K(2,3)*(grid_xP - sys.goal(3)) + sys_joint.K(2,4)*(grid_xP_dot - sys.goal(4)));
policyT_bounded = min(sys.u_limits(2,2), max(sys.u_limits(2,1), policyT));

sys_joint.V_LQR = computeValueGrid(sys_joint.S, sys_joint.l_point, grid_x, grid_x_dot, grid_xP, grid_xP_dot);
sys_joint.V_DP = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF, policyT, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
sys_joint.V_DP_bounded = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF_bounded, policyT_bounded, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);

X_DIMS = [1;2;3;4];
U_DIMS = [1;2];
% F First
disp("Cascaded F First");
numCompleted = 0;
sysFF = cell(2^length(X_DIMS) - 1, 1);
sysTS = cell(2^length(X_DIMS) - 1, 1);

for kk=1:1:4
    C = nchoosek(X_DIMS, kk);
    for ii=1:1:size(C, 1)
        sysFF_ = sys;
        sysFF_.X_DIMS_FREE = C(ii,:)';
        sysFF_.X_DIMS_FIXED = linspace(1,4,4)';
        sysFF_.X_DIMS_FIXED(sysFF_.X_DIMS_FREE) = [];
        sysFF_.U_DIMS_FREE = [1];
        sysFF_.U_DIMS_FIXED = [2];
        
        A_ = A(sysFF_.X_DIMS_FREE, sysFF_.X_DIMS_FREE);
        B_ = B(sysFF_.X_DIMS_FREE, sysFF_.U_DIMS_FREE);
        Q_ = Q(sysFF_.X_DIMS_FREE, sysFF_.X_DIMS_FREE);
        R_ = R(sysFF_.U_DIMS_FREE, sysFF_.U_DIMS_FREE);
        sysFF_.A = A_;
        sysFF_.B = B_;
        sysFF_.Q = Q_;
        sysFF_.R = R_;
        [sysFF_.K, sysFF_.S, sysFF_.e] = lqr(A_ - sysFF_.lambda_*eye(size(A_,1))/2, B_, Q_, R_, ...
                                              zeros(size(A_,1), size(B_,2)));
        K = zeros(length(sysFF_.U_DIMS_FREE), length(X_DIMS));
        K(:,sysFF_.X_DIMS_FREE) = sysFF_.K;
        policyF = -(K(1,1)*(grid_x - sys.goal(1)) + K(1,2)*(grid_x_dot - sys.goal(2)) ...
                    + K(1,3)*(grid_xP - sys.goal(3)) + K(1,4)*(grid_xP_dot - sys.goal(4)));
        policyF_bounded = min(sys.u_limits(1,2), max(sys.u_limits(1,1), policyF));
        sysFF{numCompleted + ii, 1} = sysFF_;
        
        sysTS_ = sys;
        sysTS_.X_DIMS_FREE = X_DIMS;
        sysTS_.X_DIMS_FIXED = [];
        sysTS_.U_DIMS_FREE = [2];
        sysTS_.U_DIMS_FIXED = [1];
        
        A_ = A;
        A_(:, sysFF_.X_DIMS_FREE) = A_(:, sysFF_.X_DIMS_FREE) - B(:, sysFF_.U_DIMS_FREE)*sysFF_.K;
        B_ = B;
        B_(:, sysFF_.U_DIMS_FREE) = zeros(size(B,1), length(sysFF_.U_DIMS_FREE));
        Q_ = Q;
        Q_(sysFF_.X_DIMS_FREE, sysFF_.X_DIMS_FREE) = Q_(sysFF_.X_DIMS_FREE, sysFF_.X_DIMS_FREE) ...
                                                     + sysFF_.K'*sysFF_.R*sysFF_.K;
        R_ = R;
        sysTS_.A = A_;
        sysTS_.B = B_;
        sysTS_.Q = Q_;
        sysTS_.R = R_;
        [sysTS_.K, sysTS_.S, sysTS_.e] = lqr(A_ - sysTS_.lambda_*eye(size(A_,1))/2, B_, Q_, R_, ...
                                              zeros(size(A_,1), size(B_,2)));
        sysTS_.K = sysTS_.K(sysTS_.U_DIMS_FREE, :);
        policyT = -(sysTS_.K(1,1)*(grid_x - sys.goal(1)) + sysTS_.K(1,2)*(grid_x_dot - sys.goal(2)) ...
                    + sysTS_.K(1,3)*(grid_xP - sys.goal(3)) + sysTS_.K(1,4)*(grid_xP_dot - sys.goal(4)));
        policyT_bounded = min(sys.u_limits(2,2), max(sys.u_limits(2,1), policyT));
        
        if (any(eig(sysTS_.S) < 0))
            sysTS_.S = inf*eye(size(sysTS_.S, 1));
            sysTS_.V_LQR = inf*ones(size(grid_x));
        else
            sysTS_.V_LQR = computeValueGrid(sysTS_.S, sysTS_.l_point, grid_x, grid_x_dot, grid_xP, grid_xP_dot);
        end
        
        sysTS_.V_DP = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF, policyT, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysTS_.V_DP_bounded = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                               policyF_bounded, policyT_bounded, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysTS{numCompleted + ii, 1} = sysTS_;
        disp(strcat("FF Num dims :", num2str(kk), ", Decomposition :", num2str(ii)));
        
    end
    numCompleted = numCompleted + size(C, 1);

end

% T First
disp("Cascaded T First");
numCompleted = 0;
sysTF = cell(2^length(X_DIMS) - 1, 1);
sysFS = cell(2^length(X_DIMS) - 1, 1);

for kk=1:1:4
    C = nchoosek(X_DIMS, kk);
    for ii=1:1:size(C, 1)
        sysTF_ = sys;
        sysTF_.X_DIMS_FREE = C(ii,:)';
        sysTF_.X_DIMS_FIXED = linspace(1,4,4)';
        sysTF_.X_DIMS_FIXED(sysTF_.X_DIMS_FREE) = [];
        sysTF_.U_DIMS_FREE = [2];
        sysTF_.U_DIMS_FIXED = [1];
        
        A_ = A(sysTF_.X_DIMS_FREE, sysTF_.X_DIMS_FREE);
        B_ = B(sysTF_.X_DIMS_FREE, sysTF_.U_DIMS_FREE);
        Q_ = Q(sysTF_.X_DIMS_FREE, sysTF_.X_DIMS_FREE);
        R_ = R(sysTF_.U_DIMS_FREE, sysTF_.U_DIMS_FREE);
        sysTF_.A = A_;
        sysTF_.B = B_;
        sysTF_.Q = Q_;
        sysTF_.R = R_;
        [sysTF_.K, sysTF_.S, sysTF_.e] = lqr(A_ - sysTF_.lambda_*eye(size(A_,1))/2, B_, Q_, R_, ...
                                              zeros(size(A_,1), size(B_,2)));
        K = zeros(length(sysTF_.U_DIMS_FREE), length(X_DIMS));
        K(:,sysTF_.X_DIMS_FREE) = sysTF_.K;
        policyT = -(K(1,1)*(grid_x - sys.goal(1)) + K(1,2)*(grid_x_dot - sys.goal(2)) ...
                    + K(1,3)*(grid_xP - sys.goal(3)) + K(1,4)*(grid_xP_dot - sys.goal(4)));
        policyT_bounded = min(sys.u_limits(2,2), max(sys.u_limits(2,1), policyT));
        sysTF{numCompleted + ii, 1} = sysTF_;
        
        sysFS_ = sys;
        sysFS_.X_DIMS_FREE = X_DIMS;
        sysFS_.X_DIMS_FIXED = [];
        sysFS_.U_DIMS_FREE = [1];
        sysFS_.U_DIMS_FIXED = [2];
        
        A_ = A;
        A_(:, sysTF_.X_DIMS_FREE) = A_(:, sysTF_.X_DIMS_FREE) - B(:, sysTF_.U_DIMS_FREE)*sysTF_.K;
        B_ = B;
        B_(:, sysTF_.U_DIMS_FREE) = zeros(size(B,1), length(sysTF_.U_DIMS_FREE));
        Q_ = Q;
        Q_(sysTF_.X_DIMS_FREE, sysTF_.X_DIMS_FREE) = Q_(sysTF_.X_DIMS_FREE, sysTF_.X_DIMS_FREE) ...
                                                     + sysTF_.K'*sysTF_.R*sysTF_.K;
        R_ = R;
        sysFS_.A = A_;
        sysFS_.B = B_;
        sysFS_.Q = Q_;
        sysFS_.R = R_;
        [sysFS_.K, sysFS_.S, sysFS_.e] = lqr(A_ - sysFS_.lambda_*eye(size(A_,1))/2, B_, Q_, R_, ...
                                              zeros(size(A_,1), size(B_,2)));
        sysFS_.K = sysFS_.K(sysFS_.U_DIMS_FREE, :);
        policyF = -(sysFS_.K(1,1)*(grid_x - sys.goal(1)) + sysFS_.K(1,2)*(grid_x_dot - sys.goal(2)) ...
                    + sysFS_.K(1,3)*(grid_xP - sys.goal(3)) + sysFS_.K(1,4)*(grid_xP_dot - sys.goal(4)));
        policyF_bounded = min(sys.u_limits(1,2), max(sys.u_limits(1,1), policyF));
        
        if (any(eig(sysFS_.S) < 0))
            sysFS_.S = inf*eye(size(sysFS_.S, 1));
            sysFS_.V_LQR = inf*ones(size(grid_x));
        else
            sysFS_.V_LQR = computeValueGrid(sysFS_.S, sysFS_.l_point, grid_x, grid_x_dot, grid_xP, grid_xP_dot);
        end
        
        sysFS_.V_DP = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF, policyT, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysFS_.V_DP_bounded = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                               policyF_bounded, policyT_bounded, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysFS{numCompleted + ii, 1} = sysFS_;
        disp(strcat("TF Num dims :", num2str(kk), ", Decomposition :", num2str(ii)));
    end
    numCompleted = numCompleted + size(C, 1);

end

%% Decoupled
sysDec = cell(2^length(X_DIMS) - 2, 1);
for ii=1:1:(size(sysFF,1)-1)
    sysFF_ = sysFF{ii,1};
    sysDec_ = sysFF_;
    T_DIMS_FREE = X_DIMS;
    T_DIMS_FREE(sysFF_.X_DIMS_FREE) = [];
    for jj=1:1:size(sysTF,1)
        sysTF_ = sysTF{jj,1};
        if (isempty(setdiff(T_DIMS_FREE, sysTF_.X_DIMS_FREE)) ...
            && isempty(setdiff(sysTF_.X_DIMS_FREE, T_DIMS_FREE)))
            break;
        end
    end
    
    sysDec_.K = zeros(length(U_DIMS), length(X_DIMS));
    sysDec_.K(sysFF_.U_DIMS_FREE, sysFF_.X_DIMS_FREE) = sysFF_.K;
    sysDec_.K(sysTF_.U_DIMS_FREE, sysTF_.X_DIMS_FREE) = sysTF_.K;
    
    sysDec_.A = A;
    sysDec_.B = B;
    sysDec_.Q = Q;
    sysDec_.R = R;
    sysDec_.l_point = sys.l_point;
    
    policyF = -(sysDec_.K(1,1)*(grid_x - sys.goal(1)) + sysDec_.K(1,2)*(grid_x_dot - sys.goal(2)) ...
                    + sysDec_.K(1,3)*(grid_xP - sys.goal(3)) + sysDec_.K(1,4)*(grid_xP_dot - sys.goal(4)));
    policyF_bounded = min(sys.u_limits(1,2), max(sys.u_limits(1,1), policyF));
    policyT = -(sysDec_.K(2,1)*(grid_x - sys.goal(1)) + sysDec_.K(2,2)*(grid_x_dot - sys.goal(2)) ...
                    + sysDec_.K(2,3)*(grid_xP - sys.goal(3)) + sysDec_.K(2,4)*(grid_xP_dot - sys.goal(4)));
    policyT_bounded = min(sys.u_limits(2,2), max(sys.u_limits(2,1), policyT));
    
    sysDec_.S = lyap((sysDec_.A - sysDec_.B*sysDec_.K - sysDec_.lambda_/2*eye(size(sysDec_.A,1)))'...
                      ,sysDec_.K'*sysDec_.R*sysDec_.K + sysDec_.Q);
    if (any(eig(sysDec_.S) < 0))
        sysDec_.S = inf*eye(size(sysDec_.S, 1));
        sysDec_.V_LQR = inf*ones(size(grid_x));
    else
        sysDec_.V_LQR = computeValueGrid(sysDec_.S, sysDec_.l_point, grid_x, grid_x_dot, grid_xP, grid_xP_dot);
    end
    
    sysDec_.V_DP = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF, policyT, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);                                   
    sysDec_.V_DP_bounded = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF_bounded, policyT_bounded, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
    sysDec{ii, 1} = sysDec_;
    disp(strcat("Dec, Decomposition :", num2str(ii)));
end

%% Compare
state_bounds = [-0.5, 0.5;
                -1, 1;
                2*pi/3, 4*pi/3;
                -1, 1];

valid_range = ((grid_x >= state_bounds(1,1)) & (grid_x <= state_bounds(1,2)) ...
                & (grid_x_dot >= state_bounds(2,1)) & (grid_x_dot <= state_bounds(2,2)) ...
                & (grid_xP >= state_bounds(3,1)) & (grid_xP <= state_bounds(3,2)) ...
                & (grid_xP_dot >= state_bounds(4,1)) & (grid_xP_dot <= state_bounds(4,2)));

V_decomposition_err = [];
V_err_DP = [];
V_err_DP_bounded = [];
Decompositions = cell(2*size(sysFF,1)+size(sysTF,1)-1, 2);
for ii=1:1:size(sysFF,1)
   V_decomposition_err = [V_decomposition_err; mean(abs(sysTS{ii,1}.V_LQR(valid_range) - sys_joint.V_LQR(valid_range)))];
   V_err_DP = [V_err_DP; mean(abs(sysTS{ii,1}.V_DP(valid_range) - sys_joint.V_DP(valid_range)))];
   V_err_DP_bounded = [V_err_DP_bounded; mean(abs(sysTS{ii,1}.V_DP_bounded(valid_range) - sys_joint.V_DP_bounded(valid_range)))];
   Decompositions{ii,1} = sysFF{ii,1}.X_DIMS_FREE;
   Decompositions{ii,2} = 1;
end

for ii=1:1:size(sysTF,1)
   V_decomposition_err = [V_decomposition_err; mean(abs(sysFS{ii,1}.V_LQR(valid_range) - sys_joint.V_LQR(valid_range)))];
   V_err_DP = [V_err_DP; mean(abs(sysFS{ii,1}.V_DP(valid_range) - sys_joint.V_DP(valid_range)))];
   V_err_DP_bounded = [V_err_DP_bounded; mean(abs(sysFS{ii,1}.V_DP_bounded(valid_range) - sys_joint.V_DP_bounded(valid_range)))];
   Decompositions{ii+size(sysFF,1),1} = sysTF{ii,1}.X_DIMS_FREE;
   Decompositions{ii+size(sysFF,1),2} = 2;
end

for ii=1:1:size(sysDec,1)
   V_decomposition_err = [V_decomposition_err; mean(abs(sysDec{ii,1}.V_LQR(valid_range) - sys_joint.V_LQR(valid_range)))];
   V_err_DP = [V_err_DP; mean(abs(sysDec{ii,1}.V_DP(valid_range) - sys_joint.V_DP(valid_range)))];
   V_err_DP_bounded = [V_err_DP_bounded; mean(abs(sysDec{ii,1}.V_DP_bounded(valid_range) - sys_joint.V_DP_bounded(valid_range)))];
   Decompositions{ii+size(sysFF,1)+size(sysTF,1),1} = sysDec{ii,1}.X_DIMS_FREE;
   Decompositions{ii+size(sysFF,1)+size(sysTF,1),2} = 3;
end

% V_decomposition_var = abs(V_err_DP - V_decomposition_err);
% V_decomposition_var_bounded = abs(V_err_DP_bounded - V_decomposition_err);
V_decomposition_var = (V_err_DP - V_decomposition_err);
V_decomposition_var_bounded = (V_err_DP_bounded - V_decomposition_err);

[V_decomposition_err_sorted, V_decomposition_err_order] = sort(V_decomposition_err);
for nn=1:1:size(V_decomposition_err, 1)
    decomp = Decompositions{V_decomposition_err_order(nn),1};
    decomptype = Decompositions{V_decomposition_err_order(nn),2};
    disp(strcat(num2str(decomp'), ', (', num2str(decomptype),'), Err : ', num2str(V_decomposition_err_sorted(nn)), ', Var : ', num2str(V_decomposition_var(V_decomposition_err_order(nn))), ', Var (bounded) : ', num2str(V_decomposition_var_bounded(V_decomposition_err_order(nn)))));
end

function val_grid = computeValueGrid(S, l_point, varargin)
    
    assert(size(S,1)==size(S,2), 'S must be square');
    assert(size(S,1)==(nargin-2), 'Must provide as many grid matrices as dimensions');
    assert(length(l_point)==size(S,1), 'Check l_point dimension');
    
    val_grid = zeros(size(varargin{1}));
    
    for i=1:1:size(S,1)
        for j=1:1:size(S,1)
            val_grid = val_grid + (varargin{i} - l_point(i)).*(varargin{j} - l_point(j))*S(i,j);
        end
    end
end