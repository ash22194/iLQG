function G = policyEvaluationFull(mc, mp, l, g, goal, dt, Q, R, gamma_, policyF, policyT, grid_x, grid_x_dot, grid_xP, grid_xP_dot, gtol, max_iter, varargin)
    
    if (isempty(varargin))
        G = zeros(size(grid_x));
    else
        G = varargin{1};
    end
    x_limits = [min(grid_x, [], 'all'), max(grid_x, [], 'all')];
    x_dot_limits = [min(grid_x_dot, [], 'all'), max(grid_x_dot, [], 'all')];
    xP_limits = [min(grid_xP, [], 'all'), max(grid_xP, [], 'all')];
    xP_dot_limits = [min(grid_xP_dot, [], 'all'), max(grid_xP_dot, [], 'all')];
    goal_grid = [0, 0, 0, 0];
    goal_grid(1) = find(reshape(grid_x(:,1,1,1), [size(grid_x, 1), 1])==goal(1));
    goal_grid(2) = find(reshape(grid_x_dot(1,:,1,1), [size(grid_x, 2), 1])==goal(2));
    goal_grid(3) = find(reshape(grid_xP(1,1,:,1), [size(grid_x, 3), 1])==goal(3));
    goal_grid(4) = find(reshape(grid_xP_dot(1,1,1,:), [size(grid_x, 4), 1])==goal(4));

    Cost = (Q(1,1)*(grid_x - goal(1)).^2 + Q(2,2)*(grid_x_dot - goal(2)).^2 + Q(3,3)*(grid_xP - goal(3)).^2 + Q(4,4)*(grid_xP_dot - goal(2)).^2 ...
            + (Q(1,2) + Q(2,1))*(grid_x - goal(1)).*(grid_x_dot - goal(2)) ...
            + (Q(1,3) + Q(3,1))*(grid_x - goal(1)).*(grid_xP - goal(3)) ...
            + (Q(1,4) + Q(4,1))*(grid_x - goal(1)).*(grid_xP_dot - goal(4)) ...
            + (Q(2,3) + Q(3,2))*(grid_x_dot - goal(2)).*(grid_xP - goal(3)) ...
            + (Q(2,4) + Q(4,2))*(grid_x_dot - goal(2)).*(grid_xP_dot - goal(4)) ...
            + (Q(4,3) + Q(3,4))*(grid_xP - goal(3)).*(grid_xP_dot - goal(4)))*dt ...
            + (R(1,1)*policyF.^2 + R(2,2)*policyT.^2)*dt;
    
    % -- % RK4 integration
    k1_x = dt*grid_x_dot;
    k1_x_dot = dt*(policyF - policyT.*cos(grid_xP)/l + mp*l*(grid_xP_dot.^2).*sin(grid_xP) + 0.5*mp*g*sin(2*grid_xP))./(mc + mp*sin(grid_xP).^2);
    k1_xP = dt*grid_xP_dot;
    k1_xP_dot = dt*(policyT*(mc/mp+1)/l - policyF.*cos(grid_xP) - 0.5*mp*l*(grid_xP_dot.^2).*sin(2*grid_xP) - (mc+mp)*g*sin(grid_xP))./(l*(mc + mp*sin(grid_xP).^2));
    q1_x_dot = grid_x_dot + 0.5*k1_x_dot;
    q1_xP = grid_xP + 0.5*k1_xP;
    q1_xP_dot = grid_xP_dot + 0.5*k1_xP_dot;

    k2_x = dt*q1_x_dot;
    k2_x_dot = dt*(policyF - policyT.*cos(q1_xP)/l + mp*l*(q1_xP_dot.^2).*sin(q1_xP) + 0.5*mp*g*sin(2*q1_xP))./(mc + mp*sin(q1_xP).^2);
    k2_xP = dt*q1_xP_dot;
    k2_xP_dot = dt*(policyT*(mc/mp+1)/l - policyF.*cos(q1_xP) - 0.5*mp*l*(q1_xP_dot.^2).*sin(2*q1_xP) - (mc+mp)*g*sin(q1_xP))./(l*(mc + mp*sin(q1_xP).^2));
    q2_x_dot = grid_x_dot + 0.5*k2_x_dot;
    q2_xP = grid_xP + 0.5*k2_xP;
    q2_xP_dot = grid_xP_dot + 0.5*k2_xP_dot;

    k3_x = dt*q2_x_dot;
    k3_x_dot = dt*(policyF - policyT.*cos(q2_xP)/l + mp*l*(q2_xP_dot.^2).*sin(q2_xP) + 0.5*mp*g*sin(2*q2_xP))./(mc + mp*sin(q2_xP).^2);
    k3_xP = dt*q2_xP_dot;
    k3_xP_dot = dt*(policyT*(mc/mp+1)/l - policyF.*cos(q2_xP) - 0.5*mp*l*(q2_xP_dot.^2).*sin(2*q2_xP) - (mc+mp)*g*sin(q2_xP))./(l*(mc + mp*sin(q2_xP).^2));
    q3_x_dot = grid_x_dot + k3_x_dot;
    q3_xP = grid_xP + k3_xP;
    q3_xP_dot = grid_xP_dot + k3_xP_dot;

    k4_x = dt*q3_x_dot;
    k4_x_dot = dt*(policyF - policyT.*cos(q3_xP)/l + mp*l*(q3_xP_dot.^2).*sin(q3_xP) + 0.5*mp*g*sin(2*q3_xP))./(mc + mp*sin(q3_xP).^2);
    k4_xP = dt*q3_xP_dot;
    k4_xP_dot = dt*(policyT*(mc/mp+1)/l - policyF.*cos(q3_xP) - 0.5*mp*l*(q3_xP_dot.^2).*sin(2*q3_xP) - (mc+mp)*g*sin(q3_xP))./(l*(mc + mp*sin(q3_xP).^2));

    grid_x_ = grid_x + 1/6*(k1_x + 2*k2_x + 2*k3_x + k4_x);
    grid_x_dot_ = grid_x_dot + 1/6*(k1_x_dot + 2*k2_x_dot + 2*k3_x_dot + k4_x_dot);
    grid_xP_ = grid_xP + 1/6*(k1_xP + 2*k2_xP + 2*k3_xP + k4_xP);
    grid_xP_dot_ = grid_xP_dot + 1/6*(k1_xP_dot + 2*k2_xP_dot + 2*k3_xP_dot + k4_xP_dot);
    % -- %

    greaterThan = grid_x_ > x_limits(2);
    grid_x_(greaterThan) = x_limits(2);
    lessThan = grid_x_ < x_limits(1);
    grid_x_(lessThan) = x_limits(1);

    greaterThan = grid_x_dot_ > x_dot_limits(2);
    grid_x_dot_(greaterThan) = x_dot_limits(2);
    lessThan = grid_x_dot_ < x_dot_limits(1);
    grid_x_dot_(lessThan) = x_dot_limits(1);

    % Wrap around in theta
%     greaterThan = grid_xP_ > xP_limits(2);
%     grid_xP_ = grid_xP_ - xP_limits(2)*greaterThan;
%     lessThan = grid_xP_ < xP_limits(1);
%     grid_xP_ = grid_xP_ + xP_limits(2)*lessThan;
    greaterThan = grid_xP_ > xP_limits(2);
    grid_xP_(greaterThan) = xP_limits(2);
    lessThan = grid_xP_ < xP_limits(1);
    grid_xP_(lessThan) = xP_limits(1);

    greaterThan = grid_xP_dot_ > xP_dot_limits(2);
    grid_xP_dot_(greaterThan) = xP_dot_limits(2);
    lessThan = grid_xP_dot_ < xP_dot_limits(1);
    grid_xP_dot_(lessThan) = xP_dot_limits(1);

    G_ = 10*ones(size(grid_x));
    
    iter = 0;

    while ((max(abs(G - G_),[],'all') > gtol) && (iter < max_iter))
        
        G_ = G;
        iter = iter + 1
        
        Gnext = interpn(grid_x, grid_x_dot, grid_xP, grid_xP_dot, G, grid_x_, grid_x_dot_, grid_xP_, grid_xP_dot_);
        Gnext(goal_grid(1), goal_grid(2), goal_grid(3), goal_grid(4)) = 0;
%         G = (Q(1,1)*(gridx_ - goal(1)).^2 + Q(2,2)*(gridx_dot_ - goal(2)).^2 + Q(3,3)*(gridxP_ - goal(3)).^2 + Q(4,4)*(gridxP_dot_ - goal(2)).^2 ...
%             + (Q(1,2) + Q(2,1))*(gridx_ - goal(1)).*(gridx_dot_ - goal(2)) ...
%             + (Q(1,3) + Q(3,1))*(gridx_ - goal(1)).*(gridxP_ - goal(3)) ...
%             + (Q(1,4) + Q(4,1))*(gridx_ - goal(1)).*(gridxP_dot_ - goal(4)) ...
%             + (Q(2,3) + Q(3,2))*(gridx_dot_ - goal(2)).*(gridxP_ - goal(3)) ...
%             + (Q(2,4) + Q(4,2))*(gridx_dot_ - goal(2)).*(gridxP_dot_ - goal(4)) ...
%             + (Q(4,3) + Q(3,4))*(gridxP_ - goal(3)).*(gridxP_dot_ - goal(4)))*dt ...
%             + (R(1,1)*policyF.^2 + R(2,2)*policyT.^2)*dt + gamma_*Gnext;

        G = Cost + gamma_*Gnext;
        G(goal_grid(1), goal_grid(2), goal_grid(3), goal_grid(4)) = 0;
    end
end