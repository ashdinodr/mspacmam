function flux_1d = diag_fluxes(model_grid, water_col_model, tracer_conc_1d_t, w_velocity_1d_w)
	%%% Diagnose sinking particle fluxes from tracer concentrations and w-velocities. %%%

    % Unpack select variables for convenience
    zt = model_grid.vertical.zt;
    zw = model_grid.vertical.zw;
    nzt = model_grid.vertical.nzt;
    k_max = water_col_model.vert_grid.k_max;

    % Assume that sinking fluxes are located on the edges of the vertical tracer grid cells,
    % where the w-velocities are located
    % Before the flux calculation, tracer concentrations need to be interpolated to
    % the depth coordinates on the w-velocity grid (i.e., grid cell edges)
    if k_max >= 2
        % Use linear interpolation and extrapolation methods if 1-D water column has >= 2 valid grid cells
        interp_method = 'linear'; extrap_method = 'linear';
        k_interp = k_max; % bottom-most vertical grid cell index to which interpolant should be created
    else
        % Use previous and nearest neighbor interpolation and extrapolation methods if 1-D water column has only 1 valid grid cell
        interp_method = 'previous'; extrap_method = 'nearest';
        k_interp = nzt; % N/A
    end
    interp_obj = griddedInterpolant(zt(1:k_interp),tracer_conc_1d_t(1:k_interp),interp_method,extrap_method);
    tracer_conc_1d_w = interp_obj(zw);
    % Keep only one extrapolated value
    tracer_conc_1d_w(k_max+2:end) = NaN;
    % Set any negative extrapolated values to zero
    tracer_conc_1d_w(tracer_conc_1d_w<0) = 0;

    % Sinking (advective) flux is equal to the tracer concentration times the w-velocity
    flux_1d = tracer_conc_1d_w.*w_velocity_1d_w; % mol/m^2/s
    % Display warning if any flux values are negative
    if any(flux_1d<0)
        warning('Diagnosed flux profile at (%f,%f) has one or more negative values.',water_col_model.location.lat,water_col_model.location.lon);
    end
end