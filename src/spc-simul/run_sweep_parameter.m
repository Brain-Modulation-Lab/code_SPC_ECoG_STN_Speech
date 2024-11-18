function SPC = run_sweep_parameter(cfg, param, value)
    % Run simulations by sweeping a parameter across specified values
    %
    % INPUTS:
    %   cfg    - configuration struct with simulation parameters
    %   param  - name of the parameter to sweep
    %   value  - array of values to assign to the parameter during the sweep
    %
    % OUTPUT:
    %   SPC    - struct containing the results for each value in the sweep
    
    SPC = struct(); % Initialize the output struct
    
    % Check if value is a row vector or a matrix
    if isvector(value) && (size(value, 1) == 1 || size(value, 2) == 1)
        % value is a row or column vector
        nsweeps = numel(value);
    elseif ismatrix(value) && size(value, 2) == 2
        % value is a matrix with two columns
        nsweeps = size(value, 1);
    else
        error('Value must be either a 1xN vector or an Nx2 matrix.');
    end
    
    for ii = 1:nsweeps
        % Determine the current value to use
        if isvector(value)
            currentValue = value(ii); % Sweep over vector
        else
            currentValue = value(ii, :); % Get the current row for Nx2 matrix
        end
        
        fprintf("Running simulation for param %s: ", param);
        if isvector(value)
            fprintf("%1.3f \n", currentValue);
        else
            fprintf("Row %d: [%1.3f, %1.3f] \n", ii, currentValue(1), currentValue(2));
        end
        
        cfg.(param) = currentValue;  % Set the current parameter value
        
        % Run the simulation with the updated config
        [PPC_time, PLV_time, FR, time_bin, SPC_win, cfg] = run_simulation(cfg);
        
        % Store the results for the current parameter value
        SPC(ii).param = param;
        SPC(ii).value = currentValue; % Store current value
        SPC(ii).PPC_time = PPC_time;
        SPC(ii).PLV_time = PLV_time;
        SPC(ii).FR = FR;
        SPC(ii).time_bin = time_bin;
        SPC(ii).WIN_METHOD = SPC_win;
        SPC(ii).cfg = cfg;
        
        
    end
end
