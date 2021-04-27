function number_of_decompositions = calcNumPureDecompositions(n, m)
number_of_decompositions = 0;
for r=2:1:m
%     num_r_cascaded = distr_n_into_r_with_k_nonempty(n, r, 1);
    num_r_cascaded = r^n - (r-1)^n;
    
    num_r_decoupled = distr_n_into_r(n, r);
    num_r_cascaded = num_r_cascaded * factorial(r); 
    
    num_r_groupings = distr_n_into_r(m, r);
    num_r_groupings = num_r_groupings/factorial(r);

    num_r_decompositions = num_r_groupings*(num_r_cascaded + num_r_decoupled);
    
    number_of_decompositions = number_of_decompositions + num_r_decompositions;
end

end

