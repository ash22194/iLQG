clear;
close all;
clc;

%% 

sys.m = 72;
sys.I = 3;
sys.l0 = 1.05;
sys.g = 9.81;
sys.d = 0.2;
sys.df = 0.5;
sys.dt = 0.001;
sys.X_DIMS = 6;
sys.U_DIMS = 4;
% sys.U_PSEUDO_DIMS = {[1;2];[3;4]};
sys.U_PSEUDO_DIMS = {[1];[2];[3];[4]};

lg = 0.96;
alpha1g = pi/2 + asin(sys.df/2/lg);
l2g = sqrt((sys.df + lg*cos(alpha1g))^2 + (lg*sin(alpha1g))^2);
alpha2g = acos((sys.df + lg*cos(alpha1g))/l2g);
sys.l_point = [lg; alpha1g; 0; 0; 0; 0];
sys.u0 = [sys.m*sys.g*cos(alpha2g)/sin(alpha1g - alpha2g);
          -sys.m*sys.g*cos(alpha1g)/sin(alpha1g - alpha2g);
          0;
          0];

% Continuous time dynamics
x = sym('x', [sys.X_DIMS, 1]); % l, alpha1, x_dot, z_dot, th, th_dot
u = sym('u', [sys.U_DIMS, 1]); % F_leg1, F_leg2, T_hip1, T_hip2
x_hip     = x(1)*cos(x(2));
z_hip     = x(1)*sin(x(2));
l2        = sqrt((x_hip + sys.df)^2 + z_hip^2);
a2        = acos((sys.df + x_hip)/l2);
x_hip_dot = x(3) + sys.d*cos(x(5))*x(6);
z_hip_dot = x(4) + sys.d*sin(x(5))*x(6);

zero_dyn = [0; 0; 0; -sys.g; 0; 0];
state_dyn = [x_hip_dot*cos(x(2)) + z_hip_dot*sin(x(2));
             (-x_hip_dot*sin(x(2)) + z_hip_dot*cos(x(2)))/x(1);
             0;
             0;
             x(6);
             0];
act_dyn = [0;
           0;
           (u(1)*cos(x(2)) + u(2)*cos(a2) ...
            + u(3)/x(1)*sin(x(2)) + u(4)/l2*sin(a2))/sys.m;
           (u(1)*sin(x(2)) + u(2)*sin(a2) ...
            - u(3)/x(1)*cos(x(2)) - u(4)/l2*cos(a2))/sys.m;
           0;
           (u(3)*(1 + sys.d/x(1)*sin(x(2) - x(6))) + u(4)*(1 + sys.d/l2*sin(a2 - x(6))) ...
            + u(1)*sys.d*cos(x(2) - x(6)) + u(2)*sys.d*cos(a2 - x(6)))/sys.I];
sys.f = zero_dyn + state_dyn + act_dyn;
sys.fxu = jacobian(sys.f, [x;u]);
fxu = eval(subs(sys.fxu, [x; u], [sys.l_point; sys.u0]));
sys.A = fxu(:,1:sys.X_DIMS); 
sys.B = fxu(:,(sys.X_DIMS+1):(sys.X_DIMS + sys.U_DIMS));
sys.x = x;
sys.u = u;

sys.gamma_ = 0.999;
sys.Q = diag([100, 200, 2, 2, 1000, 10]);
sys.R = 0.000002*eye(4);
sys.lambda_ = (1 - sys.gamma_)/sys.dt;

sys.numPoints = [14, 16, 17, 19, 21, 21];
sys.limits = [sys.l0 - 0.6, sys.l0 + 0.1;
              pi/2 + 0, pi/2 + 0.6;
              -0.3, 0.3;
              -0.4, 0.4;
              -pi/6, pi/6;
              -1.5, 1.5];
[grid_l, grid_a, grid_x_dot, grid_z_dot, grid_th, grid_th_dot] = ndgrid(linspace(sys.limits(1,1), sys.limits(1,2), sys.numPoints(1)), ...
                                                    linspace(sys.limits(2,1), sys.limits(2,2), sys.numPoints(2)), ...
                                                    linspace(sys.limits(3,1), sys.limits(3,2), sys.numPoints(3)), ...
                                                    linspace(sys.limits(4,1), sys.limits(4,2), sys.numPoints(4)), ...
                                                    linspace(sys.limits(5,1), sys.limits(5,2), sys.numPoints(5)), ...
                                                    linspace(sys.limits(6,1), sys.limits(6,2), sys.numPoints(6)));
sys.grid{1} = grid_l;
sys.grid{2} = grid_a;
sys.grid{3} = grid_x_dot;
sys.grid{4} = grid_z_dot;
sys.grid{5} = grid_th;
sys.grid{6} = grid_th_dot;

[K_joint, S_joint, ~] = lqr(sys.A - eye(size(sys.A,1))*sys.lambda_/2, sys.B, sys.Q, sys.R, zeros(size(sys.A,1), size(sys.B,2)));
V_joint = computeValueGrid(S_joint, sys.l_point, sys.grid{:});

state_bounds = [0.95, 1;
                pi/2 + 0.3, pi/2 + 0.4;
                -0.1, 0.1;
                -0.3, 0.3;
                -0.2, 0.2;
                -0.2, 0.2];
            
valid_range = ((grid_l >= state_bounds(1,1)) & (grid_l <= state_bounds(1,2)) ...
                & (grid_a >= state_bounds(2,1)) & (grid_a <= state_bounds(2,2)) ...
                & (grid_x_dot >= state_bounds(3,1)) & (grid_x_dot <= state_bounds(3,2)) ...
                & (grid_z_dot >= state_bounds(4,1)) & (grid_z_dot <= state_bounds(4,2)) ...
                & (grid_th >= state_bounds(5,1)) & (grid_th <= state_bounds(5,2)) ...
                & (grid_th_dot >= state_bounds(6,1)) & (grid_th_dot <= state_bounds(6,2)));

sys.valid_states = [grid_l(valid_range), ...
                    grid_a(valid_range), ...
                    grid_x_dot(valid_range), ...
                    grid_z_dot(valid_range), ...
                    grid_th(valid_range), ...
                    grid_th_dot(valid_range)]';

sys.V_joint = V_joint(valid_range);

%% Compute Measure

d = false;
r = logical([0, 0, 0, 1;
             0, 0, 0, 1;
             1, 0, 0, 0;
             1, 0, 0, 0]);
s = logical([1, 1, 1, 1, 0, 0;
             1, 1, 1, 1, 0, 0;
             0, 0, 0, 0, 1, 1;
             0, 0, 0, 0, 1, 1]);

profile on;
err_lqr_dec = computeLQRMeasureDecoupledwPseudoInputs(sys, s);
err_lqr_cas = computeLQRMeasureCascadedwPseudoInputs(sys, r, s);
profile viewer;

%% Functions

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