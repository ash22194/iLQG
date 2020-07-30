function err_lqr = computeLQRMeasureDecoupled(sys, s)
%% Inputs
% s is m x n binary matrix sij = 1 if input j is dependent on state i

s = logical(round(s));

%% Decomposition

K = zeros(sys.U_DIMS, sys.X_DIMS);
u0 = zeros(sys.U_DIMS,1);

[c, c_eq] = decoupled_constraints(s);

if (any(c_eq~=0) || any(c > 0))
    err_lqr = 1e8;
else
    for ii = 1:1:size(sys.U_PSEUDO_DIMS, 1)
    %     A_ = sys.A(s(ii,:), s(ii,:));
    %     B_ = sys.B(s(ii,:), ii);
        u0_ = u0;
        u0_(sys.U_PSEUDO_DIMS{ii}) = sys.u0(sys.U_PSEUDO_DIMS{ii});
        fxu = eval(subs(sys.fxu, [sys.x; sys.u], [sys.l_point; u0_]));
        A_ = fxu(:,1:sys.X_DIMS);
        A_ = A_(s(ii,:), s(ii,:));
        B_ = sys.B(s(ii,:), sys.U_PSEUDO_DIMS{ii}); % B calculation remains the same as system is control affine

        Q_ = sys.Q(s(ii,:), s(ii,:));
        R_ = sys.R(sys.U_PSEUDO_DIMS{ii}, sys.U_PSEUDO_DIMS{ii});

        [K(sys.U_PSEUDO_DIMS{ii}, s(ii,:)), ~, ~] = lqr(A_ - eye(size(A_,1))*sys.lambda_/2, B_, ...
                                                          Q_, R_, zeros(size(A_,1), size(B_,2)));
    end

    S = lyap((sys.A - sys.B*K - sys.lambda_/2*eye(size(sys.A,1)))', ...
              K'*sys.R*K + sys.Q);
    if (any(eig(S) < 0))
        err_lqr = inf;
    else
        %     V = computeValueGrid(S, sys.l_point, sys.grid{:});
        %     err_lqr = mean(abs(V(sys.valid_range) - sys.V_joint(sys.valid_range)), 'all');
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