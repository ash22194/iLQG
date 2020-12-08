function out = interpn_gpukernel(varargin)
% x1, ..., xn
% V
% xi1, ..., xn
% k

n = (length(varargin) - 2) / 2; % 2n + 2 inputs
assert(n == round(n), 'Check the number of inputs');

x  = varargin(1:n);
V  = varargin{n+1};
xi = varargin((n+2):(2*n+1));
k  = varargin{2*n+2}; % GPU kernel for interpolation

Nx = size(x{1}); % Number of grid points in each dimension
x1 = cellfun(@(y) y(1), x); % Min value in each dimension
dx = cellfun(@(y) y(end) - y(1), x) ./ (Nx(1:n) - 1); % Discretizatio in each dimension

% Generate corners of the n dimensional hyper-cube
corners = 0:1:(2^n - 1);
cb = zeros(n, 2^n);
for ii=1:1:n
    r = floor(corners / 2);
    cb(ii, :) = corners - 2*r;
    corners = r;
end
cb = num2cell(cb + 1, 2);
corners_index = sub2ind(Nx, cb{:}) - 1;

x1 = num2cell(x1(:));
Nx = Nx(1:n);
Nx = num2cell(Nx(:));
dx = num2cell(dx(:));

out = feval(k, xi{:}, V, Nx{:}, dx{:}, x1{:}, corners_index);
end