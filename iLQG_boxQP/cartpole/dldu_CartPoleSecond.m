function dl = dldu_CartPoleSecond(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
%     U_DIMS_FIXED = linspace(1,2,2)';
%     U_DIMS_FIXED(U_DIMS_FREE) = [];
    lims = sys.lims;
    R12 = sys.R(U_DIMS_FIXED, :);
    R12 = R12(:, U_DIMS_FREE);
    R21 = sys.R(U_DIMS_FREE, :);
    R21 = R21(:, U_DIMS_FIXED);
    R22 = sys.R(U_DIMS_FREE, :);
    R22 = R22(:, U_DIMS_FREE);
    
    dl = zeros(length(U_DIMS_FREE), size(u, 2));
    for ii=1:1:size(u, 2)
        dl(:, ii) = 2 * R22 * max(lims(U_DIMS_FREE,1),...
                                     min(lims(U_DIMS_FREE,2), ...
                                         u(:, ii))) ...
                       + (R12' + R21) * max(lims(U_DIMS_FIXED, 1), ...
                                            min(lims(U_DIMS_FIXED,2), ...
                                                k(:, ii) + K(:,:, ii)*(x(:, ii) - xn(:, ii))));
    end
    
end