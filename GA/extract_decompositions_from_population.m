function num_invalid_samples = extract_decompositions_from_population(sys, population)
    
    NUM_PSEUDO_DIMS = size(sys.U_PSEUDO_DIMS, 1);
    NUM_SAMPLES = size(population, 1);
    
    if (size(population,2) == (NUM_PSEUDO_DIMS*sys.X_DIMS))
        s = reshape((population(:, 1: NUM_PSEUDO_DIMS*sys.X_DIMS))', ...
                     NUM_PSEUDO_DIMS, sys.X_DIMS, NUM_SAMPLES);
    elseif(size(population,2) == (NUM_PSEUDO_DIMS*sys.X_DIMS + NUM_PSEUDO_DIMS^2))
        r = reshape((population(:, 1: NUM_PSEUDO_DIMS^2))', ...
                     NUM_PSEUDO_DIMS, NUM_PSEUDO_DIMS, NUM_SAMPLES);
        s = reshape((population(:, (NUM_PSEUDO_DIMS^2 + 1): (NUM_PSEUDO_DIMS^2 + NUM_PSEUDO_DIMS*sys.X_DIMS)))', ...
                     NUM_PSEUDO_DIMS, sys.X_DIMS, NUM_SAMPLES);
        
        % Check for validity
        rank_for_input = reshape(sum(r, 1), NUM_PSEUDO_DIMS, NUM_SAMPLES);
        invalid_r_samples_1 = any(rank_for_input~=1, 1);
        input_for_rank = reshape(sum(r, 2), NUM_PSEUDO_DIMS, NUM_SAMPLES);
        invalid_r_samples_2 = any(input_for_rank~=1, 1);
        invalid_r_samples = invalid_r_samples_1 | invalid_r_samples_2;
        
        state_for_inputs = reshape(sum(s, 1), sys.X_DIMS, NUM_SAMPLES);
        invalid_s_samples_1 = any(state_for_inputs~=1, 1);
        inputs_for_state = reshape(sum(s, 2), NUM_PSEUDO_DIMS, NUM_SAMPLES);
        invalid_s_samples = invalid_s_samples_1 | sum(inputs_for_state, 1)~=sys.X_DIMS;
        
        invalid_samples = invalid_r_samples | invalid_s_samples;
        
        num_invalid_samples = sum(invalid_samples);
    end
    
end