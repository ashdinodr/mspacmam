function model_grid = load_model_grid()
	%%% Read grid description from file into structure array. %%%
    
    % Load model grid information
    load('../data/input/grid_info.mat','grid');
    
    % Define horizontal grid
    model_grid.horizontal.lat = grid.lat; % latitude coordinates on horizontal tracer grid (DD)
    model_grid.horizontal.lon = grid.lon; % longitude coordinates on horizontal tracer grid (DD)
    model_grid.horizontal.nlat = length(grid.lat); % number of grid cells in latitudinal direction
    model_grid.horizontal.nlon = length(grid.lon); % number of grid cells in longitudinal direction
    % Define vertical tracer and w-velocity grids
    model_grid.vertical.zt = grid.zt; % depth coordinates on tracer grid (m)
    model_grid.vertical.zw = grid.zw; % depth coordinates on w-velocity grid (m)
    model_grid.vertical.dzt = grid.dzt; % vertical tracer grid cell thicknesses (m)
    model_grid.vertical.nzt = length(grid.zt); % number of vertical grid cell centers
    model_grid.vertical.nzw = length(grid.zw); % number of vertical grid cell faces
    % Define land-sea mask
    model_grid.landmask_3d_t = grid.lndmsk; % land-sea mask (0 for land; 1 for ocean)
    % Define bottom-most grid cell of water column (i.e., vertical grid cell directly above seafloor)
    % as well as top-most grid cell of seafloor (i.e., first vertical grid cell outside ocean model domain)
    model_grid.k_bottcell_2d_t = NaN(length(grid.lat),length(grid.lon)); % index of bottom-most ocean grid cell
    model_grid.k_topo_2d_t = NaN(length(grid.lat),length(grid.lon)); % index of top-most seafloor grid cell
    for idx_lat = 1:length(grid.lat)
        for idx_lon = 1:length(grid.lon)
            k_bottcell = find(grid.lndmsk(:,idx_lat,idx_lon)==1,1,'last');
            if k_bottcell
                model_grid.k_bottcell_2d_t(idx_lat,idx_lon) = k_bottcell;
            else
                model_grid.k_bottcell_2d_t(idx_lat,idx_lon) = 0;
            end
            k_topo = find(grid.lndmsk(:,idx_lat,idx_lon)==0,1,'first');
            if k_topo
                model_grid.k_topo_2d_t(idx_lat,idx_lon) = k_topo;
            else
                model_grid.k_topo_2d_t(idx_lat,idx_lon) = NaN;
            end
        end
    end
end