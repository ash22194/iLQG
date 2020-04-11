function dx = f_CartPoleFirst(sys, x, u)

mp = sys.mp;
mc = sys.mc;
l = sys.l;
g = sys.g;
l_point = sys.l_point;

X_DIMS_FREE = sys.X_DIMS_FREE;
X_DIMS_FIXED = sys.X_DIMS_FIXED;
% X_DIMS_FIXED = linspace(1,4,4)';
% X_DIMS_FIXED(X_DIMS_FREE) = [];

U_DIMS_FREE = sys.U_DIMS_FREE;

U = zeros(2, size(u, 2));
% lims = sys.lims(U_DIMS_FREE,:);
% U(U_DIMS_FREE, :) = max(lims(:,1)*ones(1,size(u,2)), ...
%                         min(lims(:,2)*ones(1,size(u,2)), u));
U(U_DIMS_FREE, :) = u;
X = zeros(4, size(x, 2));
X(X_DIMS_FREE, :) = x;
X(X_DIMS_FIXED, :) = l_point(X_DIMS_FIXED)*ones(1, size(x, 2));

zero_dyn = [X(2,:);
            (mp*l*(X(4,:).^2).*sin(X(3,:)) + 0.5*mp*g*sin(2*X(3,:)))./(mc + mp*sin(X(3,:)).^2);
            X(4,:);
            ((-0.5*mp*l*X(4,:).^2).*sin(2*X(3,:)) - (mc+mp)*g*sin(X(3,:)))./(l*(mc + mp*sin(X(3,:)).^2))];
act_dyn = [zeros(1, size(X,2));
           (U(1,:) - U(2,:).*cos(X(3,:))/l)./(mc + mp*sin(X(3,:)).^2);
           zeros(1, size(X,2));
           (U(2,:)/l*(mc/mp+1) - U(1,:).*cos(X(3,:)))./(l*(mc + mp*sin(X(3,:)).^2))];

dx = zero_dyn + act_dyn;
dx = dx(X_DIMS_FREE, :);

end
