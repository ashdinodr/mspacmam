function write_to_netcdf(outfile_path, var_name, dims, var_descrip, var_data)
	%%% Write specified model output variable to netCDF output file. %%%

    % Create new variable in netCDF file
    nccreate(outfile_path,var_name,'Dimensions',dims,'Format','netcdf4_classic');
    % Write data to newly created variable in netCDF file
    ncwrite(outfile_path,var_name,var_data);
    ncwriteatt(outfile_path,var_name,'Description',var_descrip);
end