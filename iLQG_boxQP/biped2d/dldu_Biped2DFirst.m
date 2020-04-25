function dl = dldu_Biped2DFirst(sys, x, u)
    
    U_DIMS_FREE = sys.U_DIMS_FREE;

    R = sys.R(U_DIMS_FREE, :);
    R = R(:, U_DIMS_FREE);
%     lims = sys.lims(U_DIMS_FREE, :);
%     dl = 2 * R * max(lims(:,1)*ones(1,size(u,2)), ...
%                      min(lims(:,2)*ones(1,size(u,2)), u));
    dl = 2 * R * u;
    
end