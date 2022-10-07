function tracer_conc_1d_w = interp_tracer(vert_grid, tracer_conc_1d_t)
	%%% Interpolate tracer concentrations to depth coordinates of grid cell edges. %%%
    
    % Interpolate tracer concentrations to depth coordinates on w-velocity grid
    if vert_grid.k_max >= 2
        % Use linear interpolation and nearest neighbor extrapolation methods if 1-D water column has >= 2 valid grid cells
        interp_method = 'linear'; extrap_method = 'nearest';
        k_interp = vert_grid.k_max; % bottom-most vertical grid cell index to which interpolant should be created
    else
        % Use previous neighbor interpolation and extrapolation methods if 1-D water column has only 1 valid grid cell
        interp_method = 'previous'; extrap_method = 'previous';
        k_interp = vert_grid.nzt; % N/A
    end
    interp_obj = griddedInterpolant(vert_grid.zt(1:k_interp),tracer_conc_1d_t(1:k_interp),interp_method,extrap_method);
    tracer_conc_1d_w = interp_obj(vert_grid.zw);
    % Keep only one extrapolated value
    tracer_conc_1d_w(vert_grid.k_max+2:end) = NaN;
    % Set any negative extrapolated values to zero
    tracer_conc_1d_w(tracer_conc_1d_w<0) = 0;
end