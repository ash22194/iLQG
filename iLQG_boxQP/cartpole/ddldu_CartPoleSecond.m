function ddl = ddldu_CartPoleSecond(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    R = sys.R(U_DIMS_FREE, :);
    R = R(:, U_DIMS_FREE);
    ddl = repmat(2 * R, [1, 1, size(x, 2)]);
end