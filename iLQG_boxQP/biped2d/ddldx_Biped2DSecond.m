function ddl = ddldx_Biped2DSecond(sys, x, u, k, K, xn)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;
    U_DIMS_FIXED = sys.U_DIMS_FIXED;
%     U_DIMS_FIXED = linspace(1,2,2)';
%     U_DIMS_FIXED(U_DIMS_FREE) = [];

    R = sys.R(U_DIMS_FIXED, :);
    R = R(:, U_DIMS_FIXED);

    ddl = zeros(6, 6, size(x, 2));
    for ii=1:1:size(x, 2)
        ddl(:,:, ii) = 2*(sys.Q + K(:,:, ii)' * R * K(:,:, ii));
    end
end