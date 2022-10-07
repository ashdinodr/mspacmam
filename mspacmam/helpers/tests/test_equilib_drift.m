function passing_result = test_equilib_drift(delta_x_1d, x_ref_1d, drift_thresh)
	%%% Test specified model variable for drift. %%%

    % First test for negative, Inf and/or NaN values
    if test_for_irregularities(x_ref_1d)
        error('Negative, Inf and/or NaN values found in steady-state solution. \nDecrease time step, dt, in config file.',[]);
    end
    
    % Test for model drift
    if all(100*(delta_x_1d./x_ref_1d) < drift_thresh)
        % To pass test, all variable values should drift by less than specified threshold
        passing_result = 1;
    else
        passing_result = 0;
    end
end