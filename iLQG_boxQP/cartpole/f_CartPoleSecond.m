function dx = f_CartPoleSecond(sys, x, u, k, K, xn)

mp = sys.mp;
mc = sys.mc;
l = sys.l;
g = sys.g;

U_DIMS_FREE = sys.U_DIMS_FREE;
U_DIMS_FIXED = sys.U_DIMS_FIXED;
%     U_DIMS_FIXED = linspace(1, 2, 2)';
%     U_DIMS_FIXED(U_DIMS_FREE) = [];

lims = sys.lims;
U = zeros(2, size(u,2));
U(U_DIMS_FREE, :) = u;
% for ii = 1:1:size(u, 2)
%     U(U_DIMS_FIXED, ii) = k(:, ii) + K(:,:, ii)*(x(:, ii) - xn(:, ii));
%     U(U_DIMS_FIXED, ii) = max(lims(U_DIMS_FIXED,1), ...
%                               min(lims(U_DIMS_FIXED,2), U(U_DIMS_FIXED, ii)));
% end

TRAJ_LENGTH = size(x,2);
x_ = x - xn;
x_ = reshape(x_, 1, size(x_,1), TRAJ_LENGTH);
x_ = repmat(x_, [length(U_DIMS_FIXED), 1, 1]);
K_ = reshape(sum(K.*x_, 2), length(U_DIMS_FIXED), TRAJ_LENGTH);
U(U_DIMS_FIXED, :) = k + K_(:,1,:);
U = max(repmat(lims(:,1), [1, size(u,2)]), ...
        min(repmat(lims(:,2), [1, size(u,2)]), U));

zero_dyn = [x(2,:);
            (mp*l*(x(4,:).^2).*sin(x(3,:)) + 0.5*mp*g*sin(2*x(3,:)))./(mc + mp*sin(x(3,:)).^2);
            x(4,:);
            ((-0.5*mp*l*x(4,:).^2).*sin(2*x(3,:)) - (mc+mp)*g*sin(x(3,:)))./(l*(mc + mp*sin(x(3,:)).^2))];
act_dyn = [zeros(1, size(x,2));
           (U(1,:) - U(2,:).*cos(x(3,:))/l)./(mc + mp*sin(x(3,:)).^2);
           zeros(1, size(x,2));
           (U(2,:)/l*(mc/mp+1) - U(1,:).*cos(x(3,:)))./(l*(mc + mp*sin(x(3,:)).^2))];

dx = zero_dyn + act_dyn;

end
