function ddl = ddldudx_Biped2DSecond(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
%     U_DIMS_FIXED = linspace(1,2,2)';
%     U_DIMS_FIXED(U_DIMS_FREE) = [];
    R12 = sys.R(U_DIMS_FIXED, :);
    R12 = R12(:, U_DIMS_FREE);
    R21 = sys.R(U_DIMS_FREE, :);
    R21 = R21(:, U_DIMS_FIXED);
    ddl = zeros(length(U_DIMS_FREE), 6, size(x, 2));
    for ii=1:1:size(x, 2)
        ddl(:,:, ii) = (K(:,:, ii)'*(R12 + R21'))';
    end
    
end