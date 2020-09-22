clear;
close all;
clc;

%% 

sys.m = 0.5;
sys.g = 9.81;
sys.l = 0.225;
sys.bk = 1.14*1e-7/(2.98*1e-6); % tau/f
sys.I = diag([4.86*1e-3; 4.86*1e-3; 8.8*1e-3]);
sys.dt = 0.001;
sys.X_DIMS = 10; % z, ro, pi, ya, vx, vy, vz, vro, vpi, vya
sys.U_DIMS = 4;
sys.U_PSEUDO_DIMS = {[1];[2];[3];[4]};

sys.l_point = [1; zeros(9,1)];
sys.u0 = sys.m*sys.g*ones(sys.U_DIMS, 1)/sys.U_DIMS;

% Continuous time dynamics
x = sym('x', [sys.X_DIMS, 1]);
u = sym('u', [sys.U_DIMS, 1]);
zero_dyn = [zeros(4,1); 0; 0; -sys.g; zeros(3, 1)];
Wn = [1, 0,          -sin(x(2));
      0, cos(x(3)),  cos(x(2))*sin(x(3));
      0, -sin(x(3)), cos(x(2))*cos(x(3))];
nu = Wn*[x(8); 
         x(9); 
         x(10)];
tau = [sys.l*(u(4) - u(2));
       sys.l*(u(3) - u(1));
       sys.bk*(u(1) - u(2) + u(3) - u(4))];
nudot = sys.I\(tau - cross(nu, sys.I*nu));

Wninvdot = [0, ...
            x(9)*cos(x(3))*tan(x(2))+x(8)*sin(x(3))/(cos(x(2))^2), ...
            -x(9)*sin(x(3))*cos(x(2))+x(8)*cos(x(3))/(cos(x(2))^2);
            0, -x(9)*sin(x(3)), -x(9)*cos(x(3));
            0, ...
            x(9)*cos(x(3))/cos(x(2))+x(9)*sin(x(3))*tan(x(2))/cos(x(2)), ...
            -x(9)*sin(x(3))/cos(x(2))+x(8)*cos(x(3))*tan(x(2))/cos(x(2))];
state_act_dyn = [x(7:10);
                 (u(1)+u(2)+u(3)+u(4))/sys.m*([cos(x(3))*sin(x(2))*cos(x(4)) + sin(x(3))*sin(x(4));
                                               cos(x(3))*sin(x(2))*sin(x(4)) - sin(x(3))*cos(x(4));
                                               cos(x(2))*cos(x(3))]);
                 Wninvdot*nu + Wn\nudot];
sys.f = zero_dyn + state_act_dyn;
sys.fxu = jacobian(sys.f, [x; u]);
sys.fxu_func = @(x, u) QuadcopterDynamicsJacobian(x, u);

fxu = eval(subs(sys.fxu, [x; u], [sys.l_point; sys.u0]));
sys.A = fxu(:,1:sys.X_DIMS); 
sys.B = fxu(:,(sys.X_DIMS+1):(sys.X_DIMS + sys.U_DIMS));
sys.x = x;
sys.u = u;
sys.xu = [sys.x; sys.u];

sys.gamma_ = 0.999;
sys.Q = diag([1, 2, 2, 0.25, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01]);
sys.R = 0.002*eye(4);
sys.lambda_ = (1 - sys.gamma_)/sys.dt;
[K_joint, S_joint, ~] = lqr(sys.A - eye(size(sys.A,1))*sys.lambda_/2, sys.B, sys.Q, sys.R, zeros(size(sys.A,1), size(sys.B,2)));
sys.S =  sym('S', [sys.X_DIMS, sys.X_DIMS]);
sys.S = tril(sys.S,0) + tril(sys.S,-1).';
sys.a = sym('a', [sys.X_DIMS, 1]);
sys.b = sym('b', [sys.X_DIMS, 1]);
sys.err_lqr = sys.x.'*sys.S*sys.x - sys.x.'*S_joint*sys.x;
for ii=1:1:sys.X_DIMS
    sys.err_lqr = int(sys.err_lqr, sys.x(ii), [sys.a(ii), -sys.a(ii)]);
end
% sys.err_lqr_func = @(S,a,b) QuadcopterErrLQR(S,a);

sys.numPoints = [12, 24*ones(1,3), 20, 20, 10, 12, 12, 10];
sys.limits = [0.4, 1.6;
              -pi*ones(3,1), pi*ones(3,1);
              -3*ones(2,1), 3*ones(2,1);
              -1.5, 1.5;
              -ones(2,1), ones(2,1);
              -0.5, 0.5];
sys.state_bounds = [0.6, 1.4;
                    -3*pi/8, 3*pi/8;
                    -3*pi/8, 3*pi/8;
                    -3*pi/8, 3*pi/8;
                    -1.5, 1.5;
                    -1.5, 1.5;
                    -0.75, 0.75;
                    -0.5, 0.5;
                    -0.5, 0.5;
                    -0.25, 0.25];
sys.da = prod(sys.state_bounds(:,2) - sys.state_bounds(:,1));

p = [0, 1;
     0, 1;
     0, 2;
     0, 2];
s = [ones(2, 2), zeros(2, 2), ones(2, 4), zeros(2, 2);
     zeros(2, 2), ones(2, 2), zeros(2, 4), ones(2, 2)];

err_lqr = computeLQRMeasure(sys, p, s);
err_compute = computeComplexityEstimates(sys, p, s);

%% Test LQR control

tspan = [0, 10];
starts = [0.4, zeros(1, 9);
          1.4, zeros(1, 9);
          1, 3*pi/8, 3*pi/8, zeros(1, 7);
          1, -3*pi/8, 3*pi/8, zeros(1, 7);
          1, zeros(1, 6), 1, 1, 0;
          1, zeros(1, 6), 1, -1, 0;
          1, zeros(1, 6), 0.5, 0.5, 0.5];
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4);
trajectories = cell(size(starts, 1), 3);

for ii=1:1:size(starts, 1)
    start = starts(ii, :)';
    [t, y] = ode45(@(t,y) QuadcopterDynamics(y, max(0, sys.u0 - K_joint*(y - sys.l_point))), tspan, start, opts);
%     [t, y] = ode45(@(t,y) QuadcopterDynamics(y, sys.u0 - K_joint*(y - sys.l_point)), tspan, start, opts);
    trajectories{ii, 1} = t;
    trajectories{ii, 2} = y;
    trajectories{ii, 3} = max(0, sys.u0 - K_joint*(y' - sys.l_point))';
    fprintf('End z:%.2f, ro:%.2f, pi:%.2f, ya:%.2f, vx:%.2f, vy:%.2f, vz:%.2f, vro:%.2f, vpi:%.2f, vya:%.2f\n', y(end, :));
end

