function k_max = find_kmax(lat_lon, model_grid, water_col_vars)
	%%% Determine index of bottom-most tracer grid cell to which model can be run. %%%

    % Extract index of bottom-most ocean grid cell at specified horizontal grid point
    horiz_grid = model_grid.horizontal;
    k_bottcell = model_grid.k_bottcell_2d_t(horiz_grid.lat==lat_lon(1),horiz_grid.lon==lat_lon(2));
    % Find index of bottom-most grid cell at specified horizontal grid point
    % below which ocean data on vertical tracer grid is incomplete
    k_notNaN = find(all(~isnan([water_col_vars.T,water_col_vars.O2,water_col_vars.omegaC,water_col_vars.omegaA]),2),1,'last');
    if isempty(k_notNaN)
        k_notNaN = 0;
    end
    % Find index of bottom-most grid cell at specified horizontal grid point
    % below which ocean data on vertical w-velocity grid is incomplete
    k_notNaN_wgrid = find(all(~isnan([water_col_vars.rho_sw_wgrid,water_col_vars.visc_wgrid]),2),1,'last')-1;
    if isempty(k_notNaN_wgrid)
        k_notNaN_wgrid = 0;
    end
    % Set k_export to zero if any export flux at specified horizontal grid point is NaN
    if any(isnan([water_col_vars.epPOC,water_col_vars.epCalc,water_col_vars.epArag,water_col_vars.epOpal]))
        k_export = 0;
    else
        k_export = Inf;
    end
    
    k_max = min([k_bottcell,k_notNaN,k_notNaN_wgrid,k_export]);
end