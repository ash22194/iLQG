function [c, c_eq] = cascaded_constraints(r, s)
    
    r = round(r);
    s = round(s);

    c = [];
    c_eq = [];
    
%     c_eq = [c_eq; (sum(r, 2) - 1)];
%     c_eq = [c_eq; (sum(r, 1)' - 1)];
%     c_eq = [c_eq; (sum(s, 1)' - 1)];
    c = [c; (1 - sum(s, 2)); 
         (sum(r, 2) - 1); (1 - sum(r, 2));
         (sum(r, 1)' - 1); (1 - sum(r, 1)');
         (sum(s, 1)' - 1); (1 - sum(s, 1)')];

end