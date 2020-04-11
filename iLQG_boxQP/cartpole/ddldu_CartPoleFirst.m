function ddl = ddldu_CartPoleFirst(sys, x, u)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    R = sys.R(U_DIMS_FREE, :);
    R = R(:, U_DIMS_FREE);
    ddl = repmat(2 * R, [1,1,size(u, 2)]);
    
end