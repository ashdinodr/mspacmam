function failing_result = test_for_irregularities(tracer_1d)
	%%% Test specified model variable for unreasonable values. %%%

    % Test for negative, Inf and/or NaN values
    if any(tracer_1d<0)
        failing_result = 1;
    elseif any(isinf(tracer_1d))
        failing_result = 1;
    elseif any(isnan(tracer_1d))
        failing_result = 1;
    else
        % To pass test, no irregular variable values should be found
        failing_result = 0;
    end
end