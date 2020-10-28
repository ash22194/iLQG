function out = dyn_e(sys, x, u)

cn = u - CN(sys, x, u);
m = M(sys, x, u);

out = zeros(8, size(x,2));
for n=1:1:size(x,2)
    out(5:8,n) = m(:,:,n) \ (cn(:, n));
end
out(1:4,:) = x(5:8,:);
end

