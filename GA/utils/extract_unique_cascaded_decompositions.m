function population_unique = extract_unique_cascaded_decompositions(sys, population)
    
    NUM_SAMPLES = size(population, 1);
    population_ = zeros(NUM_SAMPLES, sys.U_DIMS*sys.X_DIMS);
    for ii=1:1:NUM_SAMPLES
        r = population(ii, 1:sys.U_DIMS)';
        s = reshape(population(ii, (sys.U_DIMS + 1):end), sys.U_DIMS, sys.X_DIMS);
        [~, r_sorted] = sort(r);
        s_prev = [];
        for kk=1:1:sys.U_DIMS
            jj = r_sorted(kk);
            s_curr = find(s(jj,:));
            s(jj,s_prev) = 1;
            s_prev = sort(unique([s_prev, s_curr]));
        end
        population_(ii,:) = reshape(s, 1, sys.U_DIMS*sys.X_DIMS);
    end
    population_unique = unique(logical(population_), 'rows');
    population_unique = reshape(population_unique', sys.U_DIMS, sys.X_DIMS, size(population_unique,1));
end