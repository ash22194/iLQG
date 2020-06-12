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

max_iter = 10;
gtol = 0.00001;

% Joint
sys_joint = sys;
sys_joint.A = A;
sys_joint.B = B;
sys_joint.Q = Q;
sys_joint.R = R;

[sys_joint.K, sys_joint.S, sys_joint.e] = lqr(A - sys.lambda_*eye(size(A,1))/2, B, Q, R, ...
                                              zeros(size(A,1), size(B,2)));

X_DIMS = [1;2;3;4];
U_DIMS = [1;2];
% F First
disp("Cascaded F First");
numCompleted = 0;
sysFF = cell(2^length(X_DIMS) - 2, 1);
sysTS = cell(2^length(X_DIMS) - 2, 1);

for kk=1:1:3
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
        A_(:, sysFF_.X_DIMS_FREE) = A_(:, sysFF_.X_DIMS_FREE) - sysFF_.B*sysFF_.K;
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
        
        sysTS_.V_LQR = computeValueGrid(sysTS_.S, sysTS_.l_point, grid_x, grid_x_dot, grid_xP, grid_xP_dot);
        sysTS_.V_DP = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF, policyT, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysTS_.V_DP_bounded = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                               policyF_bounded, policyT_bounded, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysTS{numCompleted + ii, 1} = sysTS_;
        
    end
    numCompleted = numCompleted + size(C, 1);

end

% T First
disp("Cascaded T First");
numCompleted = 0;
sysTF = cell(2^length(X_DIMS) - 2, 1);
sysFS = cell(2^length(X_DIMS) - 2, 1);

for kk=1:1:3
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
        [K, S, sysTF_.e] = lqr(A_ - sysTF_.lambda_*eye(size(A_,1))/2, B_, Q_, R_, ...
                                              zeros(size(A_,1), size(B_,2)));
        sysTF_.K = zeros(length(sysTF_.U_DIMS_FREE), length(X_DIMS));
        sysTF_.K(:,sysTF_.X_DIMS_FREE) = K;
        sysTF_.S = zeros(length(X_DIMS));
        sysTF_.S(sysTF_.X_DIMS_FREE, sysTF_.X_DIMS_FREE) = S;
        policyT = -(sysTF_.K(1,1)*(grid_x - sys.goal(1)) + sysTF_.K(1,2)*(grid_x_dot - sys.goal(2)) ...
                    + sysTF_.K(1,3)*(grid_xP - sys.goal(3)) + sysTF_.K(1,4)*(grid_xP_dot - sys.goal(4)));
        policyT_bounded = min(sys.u_limits(1,2), max(sys.u_limits(1,1), policyT));
        sysTF{numCompleted + ii, 1} = sysTF_;
        
        sysFS_ = sys;
        sysFS_.X_DIMS_FREE = X_DIMS;
        sysFS_.X_DIMS_FIXED = [];
        sysFS_.U_DIMS_FREE = [1];
        sysFS_.U_DIMS_FIXED = [2];
        
        A_ = A;
        A_(:, sysTF_.X_DIMS_FREE) = A_(:, sysTF_.X_DIMS_FREE) - sysTF_.B*sysTF_.K;
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
        policyF_bounded = min(sys.u_limits(2,2), max(sys.u_limits(2,1), policyF));
        
        sysFS_.V_LQR = computeValueGrid(sysFS_.S, sysFS_.l_point, grid_x, grid_x_dot, grid_xP, grid_xP_dot);
        sysFS_.V_DP = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                       policyF, policyT, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysFS_.V_DP_bounded = policyEvaluationFull(sys.mc, sys.mp, sys.l, sys.g, sys.goal, sys.dt, Q, R, sys.gamma_, ...
                                               policyF_bounded, policyT_bounded, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter);
        sysFS{numCompleted + ii, 1} = sysFS_;
        disp(strcat("TF Num dims :", num2str(kk), ", Trajectory :", num2str(jj)));
    end
    numCompleted = numCompleted + size(C, 1);

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