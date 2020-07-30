clear;
close all;
clc;

%% 

sys.mc = 5; 
sys.mp = 1;
sys.l = 0.9;
sys.g = 9.81;
sys.dt = 0.001;
sys.X_DIMS = 4;
sys.U_DIMS = 2;
sys.U_PSEUDO_DIMS = {[1];[2]};

% Continuous time dynamics
x = sym('x', [sys.X_DIMS, 1]);
u = sym('u', [sys.U_DIMS, 1]);
zero_dyn = [x(2);
            (sys.mp*sys.l*(x(4)^2)*sin(x(3)) + 0.5*sys.mp*sys.g*sin(2*x(3)))/(sys.mc + sys.mp*sin(x(3))^2);
            x(4);
            (-0.5*sys.mp*sys.l*x(4)^2*sin(2*x(3)) - (sys.mc+sys.mp)*sys.g*sin(x(3)))/(sys.l*(sys.mc + sys.mp*sin(x(3))^2))];
act_dyn = [0;
           (u(1) - u(2)*cos(x(3))/sys.l)/(sys.mc + sys.mp*sin(x(3))^2);
           0;
           (u(2)/sys.l*(sys.mc/sys.mp+1) - u(1)*cos(x(3)))/(sys.l*(sys.mc + sys.mp*sin(x(3))^2))];
sys.f = zero_dyn + act_dyn;
sys.fxu = jacobian(sys.f, [x;u]);
sys.x = x;
sys.u = u;

sys.l_point = [0; 0; pi; 0];
sys.u0 = zeros(sys.U_DIMS,1);
fxu = eval(subs(sys.fxu, [x; u], [sys.l_point; sys.u0]));
sys.A = fxu(:,1:sys.X_DIMS); 
sys.B = fxu(:,(sys.X_DIMS+1):(sys.X_DIMS + sys.U_DIMS));
sys.Q = diag([25, 0.02, 25, 0.02]);
sys.R = 0.001*eye(2);
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
sys.grid{1} = grid_x;
sys.grid{2} = grid_x_dot;
sys.grid{3} = grid_xP;
sys.grid{4} = grid_xP_dot;

[K_joint, S_joint, ~] = lqr(sys.A - eye(size(sys.A,1))*sys.lambda_/2, sys.B, sys.Q, sys.R, zeros(size(sys.A,1), size(sys.B,2)));
V_joint = computeValueGrid(S_joint, sys.l_point, sys.grid{:});

state_bounds = [-0.5, 0.5;
                -1, 1;
                2*pi/3, 4*pi/3;
                -1, 1];
valid_range = ((grid_x >= state_bounds(1,1)) & (grid_x <= state_bounds(1,2)) ...
                    & (grid_x_dot >= state_bounds(2,1)) & (grid_x_dot <= state_bounds(2,2)) ...
                    & (grid_xP >= state_bounds(3,1)) & (grid_xP <= state_bounds(3,2)) ...
                    & (grid_xP_dot >= state_bounds(4,1)) & (grid_xP_dot <= state_bounds(4,2)));

sys.valid_states = [grid_x(valid_range), ...
                    grid_x_dot(valid_range), ...
                    grid_xP(valid_range), ...
                    grid_xP_dot(valid_range)]';
                
sys.V_joint = V_joint(valid_range);

%% Compute Measure

d = true;
r = logical([0, 1;
             1, 0]);
s = logical([1, 1, 0, 0;
             0, 0, 1, 1]);

err_lqr = computeLQRMeasure(sys, r, s, d);

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