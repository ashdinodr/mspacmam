function water_col_model = set_up_1d_model(lat_lon, model_grid, ocean_data, config)
	%%% Set up 1-D water column model at specified horizontal grid point. %%%
    
    % Extract 1-D water column model variables
    water_col_vars = extract_water_col(lat_lon, model_grid, ocean_data);
    	
    % Diagnose fractions of POC and calcite export, respectively,
    % that are routed to large particle size class
    sigma_poc = config.kappa_lg*max(eps,(water_col_vars.epOpal+water_col_vars.epArag)/(water_col_vars.epPOC+water_col_vars.epCalc+water_col_vars.epArag+water_col_vars.epOpal));
    sigma_calc = config.kappa_lg*max(eps,(water_col_vars.epOpal+water_col_vars.epArag)/(water_col_vars.epPOC+water_col_vars.epCalc+water_col_vars.epArag+water_col_vars.epOpal));
    
    % Determine index of bottom-most tracer grid cell to which model should be run
    k_max = find_kmax(lat_lon, model_grid, water_col_vars);
    
    water_col_model.water_col_vars = water_col_vars;
    water_col_model.water_col_vars.sigma_poc = sigma_poc;
    water_col_model.water_col_vars.sigma_calc = sigma_calc;
    
    water_col_model.vert_grid = model_grid.vertical;
    water_col_model.vert_grid.k_max = k_max;

    water_col_model.location.lat = lat_lon(1);
    water_col_model.location.lon = lat_lon(2);
end