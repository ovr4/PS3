using GLPK


function calculate_optimal_flux_distribution(stoichiometric_matrix::Array{Float64,2}, default_bounds_array::Array{Float64,2}, species_bounds_array::Array{Float64,2}, objective_coefficient_array::Array{Float64,1}; min_flag::Bool = true)

    # set some constants -
    TIME_RESTART_LIM = 60

    # Get the stoichiometric_matrix from data_dictionary -
    (number_of_species,number_of_fluxes) = size(stoichiometric_matrix);

    # # Setup the GLPK problem -
    lp_problem = GLPK.Prob();
    GLPK.set_prob_name(lp_problem, "sample");
    GLPK.set_obj_name(lp_problem, "objective")

    # Set solver parameters
    solver_parameters = GLPK.SimplexParam();
    GLPK.init_smcp(solver_parameters);
    solver_parameters.msg_lev = GLPK.MSG_OFF;

    # Are we doing min -or- max?
    if min_flag == true
    	GLPK.set_obj_dir(lp_problem, GLPK.MIN);
    else
    	GLPK.set_obj_dir(lp_problem, GLPK.MAX);
    end

    # Set the number of constraints and fluxes -
    GLPK.add_rows(lp_problem, number_of_species);
    GLPK.add_cols(lp_problem, number_of_fluxes);

    # Setup flux bounds, and objective function -
    (number_of_fluxes,number_of_bounds) = size(default_bounds_array)
    for flux_index = 1:number_of_fluxes

    	flux_lower_bound = default_bounds_array[flux_index,1]
    	flux_upper_bound = default_bounds_array[flux_index,2]

    	# Check bounds type ... default is DB -
    	if (flux_upper_bound == flux_lower_bound)
    		flux_constraint_type = GLPK.FX
    	else
    		flux_constraint_type = GLPK.DB
    	end

    	# flux symbol? (later use name - for now, fake it)
    	flux_symbol = "R_"*string(flux_index)

    	# Set the bounds in GLPK -
    	GLPK.set_col_name(lp_problem, flux_index, flux_symbol);
    	GLPK.set_col_bnds(lp_problem, flux_index, flux_constraint_type, flux_lower_bound, flux_upper_bound);
    end

    # Setup objective function -
    for (flux_index,obj_coeff) in enumerate(objective_coefficient_array)

    	# Set the objective function value in GLPK -
    	GLPK.set_obj_coef(lp_problem, flux_index, obj_coeff);
    end

    # Setup problem constraints for the metabolites -
    for species_index = 1:number_of_species

    	species_lower_bound = species_bounds_array[species_index,1]
    	species_upper_bound = species_bounds_array[species_index,2]

    	# defualt
    	species_constraint_type = GLPK.FX
    	if (species_lower_bound != species_upper_bound)
    		species_constraint_type = GLPK.DB
    	end

    	# set the symbol -
    	species_symbol = "x_"*string(species_index)

    	# Set the species bounds in GLPK -
    	GLPK.set_row_name(lp_problem, species_index, species_symbol);
    	GLPK.set_row_bnds(lp_problem, species_index, species_constraint_type, species_lower_bound, species_upper_bound);

    end

    # Setup the stoichiometric array -
    counter = 1;
    row_index_array = zeros(Int,number_of_species*number_of_fluxes);
    col_index_array = zeros(Int,number_of_species*number_of_fluxes);
    species_index_vector = collect(1:number_of_species);
    flux_index_vector = collect(1:number_of_fluxes);
    flat_stoichiometric_array = zeros(Float64,number_of_species*number_of_fluxes);
    for species_index in species_index_vector
    	for flux_index in flux_index_vector
    		row_index_array[counter] = species_index;
    		col_index_array[counter] = flux_index;
    		flat_stoichiometric_array[counter] = stoichiometric_matrix[species_index,flux_index];
    		counter+=1;
    	end
    end
    GLPK.load_matrix(lp_problem, number_of_species*number_of_fluxes, row_index_array, col_index_array, flat_stoichiometric_array);

    # Call the solver -
    exit_flag = GLPK.simplex(lp_problem, solver_parameters);

    # Get the objective function value -
    objective_value = GLPK.get_obj_val(lp_problem);

    # Get the calculated flux values from GLPK -
    calculated_flux_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	calculated_flux_array[flux_index] = GLPK.get_col_prim(lp_problem, flux_index);
    end

    # Get the dual values -
    dual_value_array = zeros(Float64,number_of_fluxes);
    for flux_index in flux_index_vector
    	dual_value_array[flux_index] = GLPK.get_col_dual(lp_problem, flux_index);
    end

    # is this solution optimal?
    status_flag = GLPK.get_status(lp_problem)

    # Calculate the uptake array -
    uptake_array = stoichiometric_matrix*calculated_flux_array;

    # Formulate the return tuple -
    return (objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag);
end
