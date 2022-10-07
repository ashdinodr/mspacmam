function output_1d = run_1d_model(water_col_model, config, constants)
	%%% Run MSPACMAM in 1-D water column at specified horizontal grid point. %%%
    
    % Unpack select variables for convenience
    dt = config.dt;
    nzt = water_col_model.vert_grid.nzt;
    k_max = water_col_model.vert_grid.k_max;
    
    % Determine number of time steps from user specified step size
    nt = ceil(config.runlen/dt);

    % Initialize 1-D arrays of shape (z,), filled with NaNs,
    % to store particle tracer concentrations
    prog_vars_1d.POCsm = NaN(nzt,1); % small particle POC concentration [POC] (mol C/m^3)
    prog_vars_1d.POClg = NaN(nzt,1); % large particle POC concentration [POC] (mol C/m^3)
    prog_vars_1d.CaCO3calcsm = NaN(nzt,1); % small particle calcite concentration [CaCO3] (mol C/m^3)
    prog_vars_1d.CaCO3calclg = NaN(nzt,1); % large particle calcite concentration [CaCO3] (mol C/m^3)
    prog_vars_1d.CaCO3aragsm = NaN(nzt,1); % small particle aragonite concentration [CaCO3] (mol C/m^3)
    prog_vars_1d.CaCO3araglg = NaN(nzt,1); % large particle aragonite concentration [CaCO3] (mol C/m^3)
    prog_vars_1d.SiO2sm = NaN(nzt,1); % small particle opal concentration [SiO2] (mol Si/m^3)
    prog_vars_1d.SiO2lg = NaN(nzt,1); % large particle opal concentration [SiO2] (mol Si/m^3)
    % Set initial value y0 to zero at initial time t0, i.e., y0 = y(t0) = 0
    prog_vars_1d.POCsm(1:k_max) = 0;
    prog_vars_1d.POClg(1:k_max) = 0;
    prog_vars_1d.CaCO3calcsm(1:k_max) = 0;
    prog_vars_1d.CaCO3calclg(1:k_max) = 0;
    prog_vars_1d.CaCO3aragsm(1:k_max) = 0;
    prog_vars_1d.CaCO3araglg(1:k_max) = 0;
    prog_vars_1d.SiO2sm(1:k_max) = 0;
    prog_vars_1d.SiO2lg(1:k_max) = 0;

    % Loop over time steps
    for idx_t = 1:nt
        if config.integ_method == 1
            % Approximate solutions to differential equations at the next time step
            % using forward Euler method, i.e., using the derivative information
            % from the beginning of the time step:
            % y_(n+1) = y_n + dt*y'(t_n), where y'(t_n) = dy/dt = f(t_n, y_n)
            
            % Evaluate the derivative function at the starting point of the time step
            deriv_func_output = eval_deriv_func(prog_vars_1d, water_col_model, config, constants);
            dCdt = deriv_func_output.dCdt; % derivative at beginning of interval

        elseif config.integ_method == 2
            % Approximate solutions to differential equations at the next time step
            % using RK2 midpoint method, i.e., using the derivative information
            % from the midpoint of the time step:
            % y_(n+1) = y_n + dt*f(t_n + dt/2, y_n + dt/2*f(t_n, y_n)),
            % where function f gives the derivative
            % RK2 algorithm can be written as:
            % k1 = f(t_n, y_n)
            % y1 = y_n + dt/2*k1
            % k2 = f(t_n + dt/2, y1)
            % y_(n+1) = y_n + dt*k2
            
            % Evaluate the derivative function at the starting point of the time step
            deriv_func_output = eval_deriv_func(prog_vars_1d, water_col_model, config, constants);
            k1 = deriv_func_output.dCdt; % derivative k1 at beginning of interval
            % Estimate y value at the midpoint using the derivative information
            [y1.POCsm, y1.POClg, y1.CaCO3calcsm, y1.CaCO3calclg, y1.CaCO3aragsm, y1.CaCO3araglg, y1.SiO2sm, y1.SiO2lg]...
                = deal(NaN(nzt,1), NaN(nzt,1), NaN(nzt,1), NaN(nzt,1), NaN(nzt,1), NaN(nzt,1), NaN(nzt,1), NaN(nzt,1));
            y1.POCsm(1:k_max) = prog_vars_1d.POCsm(1:k_max) + k1.SMS_POCsm(1:k_max)*dt/2;
            y1.POClg(1:k_max) = prog_vars_1d.POClg(1:k_max) + k1.SMS_POClg(1:k_max)*dt/2;
            y1.CaCO3calcsm(1:k_max) = prog_vars_1d.CaCO3calcsm(1:k_max) + k1.SMS_CaCO3calcsm(1:k_max)*dt/2;
            y1.CaCO3calclg(1:k_max) = prog_vars_1d.CaCO3calclg(1:k_max) + k1.SMS_CaCO3calclg(1:k_max)*dt/2;
            y1.CaCO3aragsm(1:k_max) = prog_vars_1d.CaCO3aragsm(1:k_max) + k1.SMS_CaCO3aragsm(1:k_max)*dt/2;
            y1.CaCO3araglg(1:k_max) = prog_vars_1d.CaCO3araglg(1:k_max) + k1.SMS_CaCO3araglg(1:k_max)*dt/2;
            y1.SiO2sm(1:k_max) = prog_vars_1d.SiO2sm(1:k_max) + k1.SMS_SiO2sm(1:k_max)*dt/2;
            y1.SiO2lg(1:k_max) = prog_vars_1d.SiO2lg(1:k_max) + k1.SMS_SiO2lg(1:k_max)*dt/2;

            % Evaluate the derivative function at the midpoint of the time step
            deriv_func_output = eval_deriv_func(y1, water_col_model, config, constants);
            dCdt = deriv_func_output.dCdt; % derivative k2 at midpoint of interval
        end
        
        if idx_t == nt
            % Throw error if equilibrium drift criteria are not met
            if ~test_equilib_drift(dCdt.SMS_POCsm(1:k_max)*dt, prog_vars_1d.POCsm(1:k_max), config.drift_thresh)
                error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
            end
            if ~test_equilib_drift(dCdt.SMS_POClg(1:k_max)*dt, prog_vars_1d.POClg(1:k_max), config.drift_thresh)
                error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
            end
            if ~test_equilib_drift(dCdt.SMS_CaCO3calcsm(1:k_max)*dt, prog_vars_1d.CaCO3calcsm(1:k_max), config.drift_thresh)
                error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
            end
            if ~test_equilib_drift(dCdt.SMS_CaCO3calclg(1:k_max)*dt, prog_vars_1d.CaCO3calclg(1:k_max), config.drift_thresh)
                error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
            end
            if config.fragment_flag == 1
                if ~test_equilib_drift(dCdt.SMS_CaCO3aragsm(1:k_max)*dt, prog_vars_1d.CaCO3aragsm(1:k_max), config.drift_thresh)
                    error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
                end
            end
            if ~test_equilib_drift(dCdt.SMS_CaCO3araglg(1:k_max)*dt, prog_vars_1d.CaCO3araglg(1:k_max), config.drift_thresh)
                error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
            end
            if config.fragment_flag == 1
                if ~test_equilib_drift(dCdt.SMS_SiO2sm(1:k_max)*dt, prog_vars_1d.SiO2sm(1:k_max), config.drift_thresh)
                    error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
                end
            end
            if ~test_equilib_drift(dCdt.SMS_SiO2lg(1:k_max)*dt, prog_vars_1d.SiO2lg(1:k_max), config.drift_thresh)
                error('Model drift exceeds threshold of %f%% per time step. \nIncrease integration time, runlen, in config file.',config.drift_thresh);
            end
        end

        % Determine the next value y_(n+1) using the derivative information
        prog_vars_1d.POCsm(1:k_max) = prog_vars_1d.POCsm(1:k_max) + dCdt.SMS_POCsm(1:k_max)*dt;
        prog_vars_1d.POClg(1:k_max) = prog_vars_1d.POClg(1:k_max) + dCdt.SMS_POClg(1:k_max)*dt;
        prog_vars_1d.CaCO3calcsm(1:k_max) = prog_vars_1d.CaCO3calcsm(1:k_max) + dCdt.SMS_CaCO3calcsm(1:k_max)*dt;
        prog_vars_1d.CaCO3calclg(1:k_max) = prog_vars_1d.CaCO3calclg(1:k_max) + dCdt.SMS_CaCO3calclg(1:k_max)*dt;
        prog_vars_1d.CaCO3aragsm(1:k_max) = prog_vars_1d.CaCO3aragsm(1:k_max) + dCdt.SMS_CaCO3aragsm(1:k_max)*dt;
        prog_vars_1d.CaCO3araglg(1:k_max) = prog_vars_1d.CaCO3araglg(1:k_max) + dCdt.SMS_CaCO3araglg(1:k_max)*dt;
        prog_vars_1d.SiO2sm(1:k_max) = prog_vars_1d.SiO2sm(1:k_max) + dCdt.SMS_SiO2sm(1:k_max)*dt;
        prog_vars_1d.SiO2lg(1:k_max) = prog_vars_1d.SiO2lg(1:k_max) + dCdt.SMS_SiO2lg(1:k_max)*dt;
    end
    
    % Throw error if there are any irregularities in the steady-state solutions
    if test_for_irregularities(prog_vars_1d.POCsm(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    if test_for_irregularities(prog_vars_1d.POClg(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    if test_for_irregularities(prog_vars_1d.CaCO3calcsm(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    if test_for_irregularities(prog_vars_1d.CaCO3calclg(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    if test_for_irregularities(prog_vars_1d.CaCO3aragsm(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    if test_for_irregularities(prog_vars_1d.CaCO3araglg(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    if test_for_irregularities(prog_vars_1d.SiO2sm(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    if test_for_irregularities(prog_vars_1d.SiO2lg(1:k_max))
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    
    output_1d.prog_vars_1d = prog_vars_1d;
    output_1d.diag_vars_1d = deriv_func_output.diag_vars_1d;
end