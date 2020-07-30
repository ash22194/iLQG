function err_lqr = computeLQRMeasureCascadedwRealRanks(sys, r, s)
%% Inputs
% r is m vector where rj is the rank for input j
% s is m x n binary matrix sij = 1 if input j is dependent on state i

r = round(r);
s = logical(round(s));

%% Decomposition
K = zeros(sys.U_DIMS, sys.X_DIMS);
u0 = zeros(sys.U_DIMS, 1);

% Check constraints
[c, c_eq] = cascaded_constraints_wreal_ranks(r, s);

if (any(c_eq~=0) || any(c > 0))
    err_lqr = 1e8;
else
    U_PSEUDO_DIMS = {};
    for jj=1:1:sys.U_DIMS
        u_r = find(r==jj);
        if (~isempty(u_r))
            U_PSEUDO_DIMS{end+1,1} = u_r;
        end
    end

    u_prev = [];
    s_curr = [];

    for ii=1:1:size(U_PSEUDO_DIMS, 1)
        u_curr = U_PSEUDO_DIMS{ii,1};
        u_prev = sort([u_prev; u_curr]);
        s_curr = sort([s_curr, find(s(u_curr(1),:))]);

    %         A_ = sys.A - sys.B*K;
        u0_ = u0;
        u0_(u_prev) = sys.u0(u_prev);
        fxu = eval(subs(sys.fxu, [sys.x; sys.u], [sys.l_point; u0_]));
        A_ = fxu(:,1:sys.X_DIMS);
        A_ = A_ - sys.B*K; % B calculation remains the same as system is control affine
        A_ = A_(s_curr, s_curr);
        B_ = zeros(length(s_curr), size(sys.B,2)); 
        B_(:, u_curr) = sys.B(s_curr, u_curr);
        R_ = sys.R;
        Q_ = sys.Q + K'*R_*K;
        
        try
            [K_, S, ~] = lqr(A_ - eye(size(A_,1))*sys.lambda_/2, B_, ...
                             Q_, R_, zeros(size(A_,1), size(B_,2)));
            K(u_curr, s_curr) = K_(u_curr, :);
        catch ME
            S = -ones(sys.X_DIMS);
            break;
%             disp(ME.identifier);
        end
        
    end

    if (any(eig(S)< 0))
        err_lqr = inf;
    else
    %     V = computeValueGrid(S, sys.l_point, sys.grid{:});
    %     err_lqr = mean(abs(V(sys.valid_range) - sys.V_joint(sys.valid_range)), 'all');
        if (size(S)~=6)
            disp('something is wrong!');
        end
        V = sum((sys.valid_states - sys.l_point).*(S*(sys.valid_states - sys.l_point)), 1)';
        err_lqr = mean(abs(V - sys.V_joint));
    end
end
end

function val_grid = computeValueGrid(S, l_point, varargin)
    
    assert(size(S,1)==size(S,2), 'S must be square');
    assert(size(S,1)==(nargin-2), 'Must provide as many grid matrices as dimensions');
    assert(length(l_point)==size(S,1), 'Check l_point dimension');
    
    val_grid = zeros(size(varargin{1}));
    
    for i=1:1:size(S,1)
        for j=1:1:size(S,1)
            val_grid = val_grid + (varargin{i} - l_point(i)).*(varargin{j} - l_point(j))*S(i,j);
        end
    end
end