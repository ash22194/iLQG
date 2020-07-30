function [c, c_eq] = decoupled_constraints_wpseudo_inputs(s)
    
    s = round(s);
    m = size(s, 1);
    
    c = [];
    c_eq = [];
    
    c = [c; (1 - sum(s, 2)); (1 - sum(s,1)')];
    
    input_couples = nchoosek(linspace(1,m,m), 2);
    for ii=1:1:size(input_couples,1)
        complete_state_overlap = 1 - min(1, sum(abs(s(input_couples(ii,1),:) - s(input_couples(ii,2),:))));
        no_state_overlap = 1 - min(1, sum(s(input_couples(ii,1),:).*s(input_couples(ii,2),:)));
        
%         c_eq = [c_eq; no_state_overlap*(1 - complete_state_overlap) ...
%                       + complete_state_overlap*(1 - no_state_overlap) - 1];
          c = [c; no_state_overlap*(1 - complete_state_overlap) ...
                  + complete_state_overlap*(1 - no_state_overlap) - 1;
                  1 - no_state_overlap*(1 - complete_state_overlap) ...
                  - complete_state_overlap*(1 - no_state_overlap)];
    end
end