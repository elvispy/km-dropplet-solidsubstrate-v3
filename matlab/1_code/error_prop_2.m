
function sigma_f_func = error_prop_2(expr, variables)
    % Create symbolic variables for uncertainties
    uncertainties = sym("s" + string(variables));
    
    % Compute the uncertainty propagation formula
    terms = arrayfun(@(var, svar) diff(expr, var)^2 * svar^2, ...
                     variables, uncertainties, 'UniformOutput', false);
    sigma_f = sqrt(sum(cat(1, terms{:}))); % Concatenate cell contents into an array
    % Combine variables and their uncertainties
    all_vars = [variables, uncertainties];
    
    % Generate a MATLAB function using matlabFunction
    sigma_f_func = matlabFunction(sigma_f, 'Vars', all_vars);
end