function [c, c_eq] = decoupled_constraints(s)
    
    s = round(s);
    
    c = [];
    c_eq = [];
    
    c = [c; (1 - sum(s, 2));
         (sum(s, 1)' - 1); (1 - sum(s, 1)')];
%     c_eq = [c_eq; (sum(s, 1)' - 1)];
    
end