function [PLV_across_trials, PPC_across_trials, nspikes] = calculate_PLV_traditional(spike_train_all_trials, osc_signal_all_trials, time_windows, t_trial, n_trials)
    % Preallocate storage for results
    n_windows = size(time_windows, 1);                  % Number of time windows
    PLV_across_trials = zeros(1, n_windows);            % Initialize PLV for each time window
    PPC_across_trials = zeros(1, n_windows);            % Initialize PPC for each time window
    nspikes = nan(1, n_windows);                         % Initialize number of spikes for each window

    % Loop over each time window
    for w = 1:n_windows
        % Define the current time window
        window_start = mean(time_windows(w, :)) - 0.1;  % Start time of the window
        window_end = mean(time_windows(w, :)) + 0.1;    % End time of the window

        % Find indices in the trial time vector that correspond to the current window
        window_idx = (t_trial >= window_start) & (t_trial < window_end);

        % Initialize array to accumulate spike phases across all trials
        spike_phases_across_trials = [];

        % Loop over each trial to collect spike phases
        for trial = 1:n_trials
            spike_train = spike_train_all_trials(trial, :);                % Get spike train for current trial
            osc_phase = angle(hilbert(osc_signal_all_trials(trial, :)));   % Compute phase of the oscillatory signal

            % Get spike times and phases within the current window using logical indexing
            spikes_in_window = spike_train(window_idx);  % Logical array of spikes in this window
            phases_in_window = osc_phase(window_idx);     % Phases corresponding to the time window

            % Collect phases of spikes within this time window
            spike_phases_across_trials = [spike_phases_across_trials, phases_in_window(spikes_in_window == 1)];
        end

        % Calculate Phase Locking Value (PLV) for the current window
        if ~isempty(spike_phases_across_trials)
            PLV_across_trials(w) = abs(mean(exp(1i * spike_phases_across_trials))); % Calculate the PLV
        else
            PLV_across_trials(w) = 0; % If no spikes in the window, set PLV to 0
        end

        % Count the number of spikes that occurred in this window
        nspikes(w) = numel(spike_phases_across_trials);

        % Calculate Phase Locking Value adjusted for number of spikes (PPC)
        PPC_across_trials(w) = correct_PLV(PLV_across_trials(w), nspikes(w)); % Adjusted PLV

        if w == 37
            figure
            histogram(spike_phases_across_trials)
        end

    end

    % Set any negative PPC values to zero
    PPC_across_trials(PPC_across_trials < 0) = 0; 
end
