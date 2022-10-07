function config = read_config_file(constants)
	%%% Read config file into structure array. %%%
	
    % Run 'config.m' file
    run('config.m');
    % List workspace variables
    config_vars = who();
    
    % Convert units if needed
    dt = round(dt*constants.secInDay); % d -> s
    runlen = runlen*constants.secInDay; % d -> s
    r_sm = r_sm/10000; % micron -> cm
    r_lg = r_lg/10000; % micron -> cm
    k_POC = k_POC/constants.secInDay; % 1/d -> 1/s
    aE = log(Q10)/10; % 1/degC
    K_O2 = K_O2/1000000/0.001; % micromol/L -> mol/m^3
    
    config = struct();
    for k=1:length(config_vars)
        config = setfield(config,config_vars{k},eval(config_vars{k}));
    end
    config = setfield(config,'aE',aE);
    config = rmfield(config,'Q10');
    config = rmfield(config,'constants');
end