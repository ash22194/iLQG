function fu = fu_CartPoleFirst(sys, x, u)

mp = sys.mp;
mc = sys.mc;
l = sys.l;
l_point = sys.l_point;
X_DIMS_FREE = sys.X_DIMS_FREE;
X_DIMS_FIXED = sys.X_DIMS_FIXED;
% X_DIMS_FIXED = linspace(1,4,4)';
% X_DIMS_FIXED(X_DIMS_FREE) = [];

U_DIMS_FREE = sys.U_DIMS_FREE;

X = zeros(4, size(x, 2));
X(X_DIMS_FREE, :) = x;
X(X_DIMS_FIXED, :) = l_point(X_DIMS_FIXED)*ones(1, size(x, 2));

fu = zeros(4, 2, size(x, 2));
fu(2,1,:) = 1.0./(mc + mp*sin(X(3,:)).^2);
fu(2,2,:) = -cos(X(3,:))./(l*(mc + mp*sin(X(3,:)).^2));
fu(4,1,:) = -cos(X(3,:))./(l*(mc + mp*sin(X(3,:)).^2));
fu(4,2,:) = (1/l*(mc/mp+1))./(l*(mc + mp*sin(X(3,:)).^2));

fu = fu(X_DIMS_FREE, :, :);
fu = fu(:, U_DIMS_FREE, :);
end
