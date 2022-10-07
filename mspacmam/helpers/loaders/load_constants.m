function constants = load_constants()
	%%% Read constants from file into structure array. %%%
	
    % Run 'constants.m' file
    run('constants.m');
    % List workspace variables
    constant_vars = who();
    
    % Convert units if needed
    g = g*100; % m/s^2 -> cm/s^2
    k_calcUp = k_calcUp*10000; % mol/cm^2/s -> mol/m^2/s
    k_calcLow = k_calcLow*10000; % mol/cm^2/s -> mol/m^2/s
    k_arag = k_arag/2.46/100/secInDay; % g/g/d -> mol/m^2/s
    
    constants = struct();
    for k=1:length(constant_vars)
        constants = setfield(constants,constant_vars{k},eval(constant_vars{k}));
    end
end