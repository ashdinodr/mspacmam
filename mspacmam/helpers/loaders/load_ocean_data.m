function ocean_data = load_ocean_data()
	%%% Read ocean state variables from file into structure array. %%%

    % Load prescribed ocean state variables
    load('../data/input/ocean_data.mat','driving_vars');
    % Load prescribed particle export fluxes
    load('../data/input/export_data.mat','export');
    
    % Convert units if needed
    driving_vars.rho_w = driving_vars.rho_w/1000; % kg/m^3 -> g/cm^3
    driving_vars.visc_w = driving_vars.visc_w*10; % kg/m/s -> g/cm/s
    
    ocean_data.driving_vars = driving_vars;
    ocean_data.export = export;
end