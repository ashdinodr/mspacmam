function valid_lats_lons = load_user_gridpoints()
	%%% Read horizontal grid points from user given file into 2-D array. %%%
    
    % Load user specified horizontal grid points
    fileID = fopen('my_ocean_gridpoints.txt','r');
    textscan(fileID,'%f %f','CommentStyle','//');
    C_coords = textscan(fileID,'%f %f','CollectOutput',1,'HeaderLines',1,'Delimiter','\t');
    fclose(fileID);
    user_lats_lons = C_coords{:};

    % Bin grid points to model grid
    load('../data/input/grid_info.mat','grid');
    [~,~,~,bin_lat,bin_lon] = histcounts2(user_lats_lons(:,1),user_lats_lons(:,2),grid.lat_b,grid.lon_b);
    rebinned_lats_lons = NaN(size(user_lats_lons));
    for idx_coord = 1:size(rebinned_lats_lons,1)
        rebinned_lats_lons(idx_coord,:) = [mean(grid.lat_b(bin_lat(idx_coord):bin_lat(idx_coord)+1)),...
            mean(grid.lon_b(bin_lon(idx_coord):bin_lon(idx_coord)+1))];
    end

    % Remove grid points outside ocean model domain
    load('../data/input/all_ocean_gridpoints.mat','lat_lon_coords');
    valid_lats_lons = [];
    for idx_coord = 1:size(rebinned_lats_lons,1)
        if ismember(rebinned_lats_lons(idx_coord,:),lat_lon_coords,'rows')
            valid_lats_lons(end+1,:) = rebinned_lats_lons(idx_coord,:);
        end
    end
end