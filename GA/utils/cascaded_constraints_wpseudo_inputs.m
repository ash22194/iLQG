function [c, c_eq] = cascaded_constraints_wpseudo_inputs(r, s)
    
    r = round(r);
    s = round(s);

    m = size(s,1);
    
    c = [];
    c_eq = [];
    
%     c_eq = [c_eq; (sum(r, 2) - 1)];
    c = [c; (1 - sum(s, 2)); (1 - sum(s, 1)');
         (sum(r, 2) - 1); (1 - sum(r, 2))];
    
    input_couples = nchoosek(linspace(1,m,m), 2);
    for ii=1:1:size(input_couples,1)
        are_coupled = min(1, sum(prod(r(input_couples(ii,:),:),1)));
        complete_state_overlap = min(1, sum(abs(s(input_couples(ii,1),:) - s(input_couples(ii,2),:))));
        no_state_overlap = min(1, sum(s(input_couples(ii,1),:).*s(input_couples(ii,2),:)));
        
%         c_eq = [c_eq; are_coupled*complete_state_overlap];
%         c_eq = [c_eq; (1-are_coupled)*no_state_overlap];
        c = [c; are_coupled*complete_state_overlap; -are_coupled*complete_state_overlap;
             (1-are_coupled)*no_state_overlap; (are_coupled-1)*no_state_overlap];
    end
    
    % Constraint to avoid jointly optimizing the inputs
    c = [c; sum(s, 2) - (sys.X_DIMS-1)];
end