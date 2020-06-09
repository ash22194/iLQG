function [xnew,unew,cnew,kFclose,KFclose,xFclose] = forward_pass_second_kdtree(x0,u,L,x,du,kF,KF,xF,Mdl,X_DIMS_FIRST,Alpha,gamma_,DYNCST,lims,diff)
% parallel forward-pass (rollout)
% internally time is on the 3rd dimension, 
% to facillitate vectorized dynamics calls

n        = size(x0,1);
K        = length(Alpha);
K1       = ones(1,K); % useful for expansion
m        = size(u,1);
N        = size(u,2);

xnew        = zeros(n,K,N);
xnew(:,:,1) = x0(:,ones(1,K));
unew        = zeros(m,K,N);
cnew        = zeros(1,K,N+1);
kFclose     = zeros(size(kF, 1), K, N+1);
KFclose     = zeros(size(KF, 1), size(KF, 2), K, N+1);
xFclose     = zeros(n,K,N+1);

discount = 1;
for i = 1:N
    unew(:,:,i) = u(:,i*K1);
    
    if ~isempty(du)
        unew(:,:,i) = unew(:,:,i) + du(:,i)*Alpha;
    end    
    
    if ~isempty(L)
        if ~isempty(diff)
            dx = diff(xnew(:,:,i), x(:,i*K1));
        else
            dx          = xnew(:,:,i) - x(:,i*K1);
        end
        unew(:,:,i) = unew(:,:,i) + L(:,:,i)*dx;
    end
    
    if ~isempty(lims)
        unew(:,:,i) = min(lims(:,2*K1), max(lims(:,1*K1), unew(:,:,i)));
    end
%     closest_xn2 = knnsearch(Mdl, xnew(X_DIMS_FIRST,:,i)');
%     xFclose(X_DIMS_FIRST,:,i) = Mdl.X(closest_xn2, :)';
    [ closest_xn2, ~, Mdl ] = kdtreeidx([], xnew(X_DIMS_FIRST,:,i)', Mdl);
    xFclose(X_DIMS_FIRST,:,i) = xF(:, closest_xn2);
    kFclose(:,:,i) = kF(:, closest_xn2);
    KFclose(:,:,:,i) = KF(:,:,closest_xn2);
    
    [xnew(:,:,i+1), cnew(:,:,i)]  = DYNCST(xnew(:,:,i), unew(:,:,i), ...
                                           kFclose(:,:,i), KFclose(:,:,:,i), xFclose(:,:,i), ...
                                           i*K1);
    cnew(:,:,i) = discount*cnew(:,:,i);
    discount = discount*gamma_;
end
% closest_xn2 = knnsearch(Mdl, xnew(X_DIMS_FIRST,:,N+1)');
% xFclose(X_DIMS_FIRST,:,N+1) = Mdl.X(closest_xn2,:)';
[ closest_xn2, ~, Mdl ] = kdtreeidx([], xnew(X_DIMS_FIRST,:,N+1)', Mdl);
xFclose(X_DIMS_FIRST,:,N+1) = xF(:, closest_xn2);
kFclose(:,:,N+1) = kF(:, closest_xn2);
KFclose(:,:,:,N+1) = KF(:,:,closest_xn2);

[~, cnew(:,:,N+1)] = DYNCST(xnew(:,:,N+1),nan(m,K,1), ...
                            kFclose(:,:,N+1), KFclose(:,:,:,N+1), xFclose(:,:,N+1), ...
                            i);
cnew(:,:,N+1) = discount*cnew(:,:,N+1);
% put the time dimension in the columns
xnew = permute(xnew, [1 3 2]);
unew = permute(unew, [1 3 2]);
cnew = permute(cnew, [1 3 2]);
xFclose = permute(xFclose, [1 3 2]);
kFclose = permute(kFclose, [1 3 2]);
KFclose = permute(KFclose, [1 2 4 3]);
end % forward pass
