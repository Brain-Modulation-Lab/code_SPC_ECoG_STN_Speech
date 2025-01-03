function [PPC_across_trials, PLV_across_trials, FR, time_bin, SPC_method,cfg] = run_simulation(cfg)
    % Run simulations for phase coupling based on given configuration.
    %
    % INPUTS:
    %   cfg             - configuration struct with simulation parameters
    %
    % OUTPUTS:
    %   PPC_across_trials - Phase-Peak Coupling results across trials
    %   PLV_across_trials - Phase Locking Value results across trials
    %   FR               - Firing rate across trials
    %   time_bin         - Time bins for results
    %   SPC_method       - Spike-phase coupling method results
    
    % Unpack configuration parameters
    fs = cfg.fs;
    trial_duration = cfg.trial_duration; % Trial duration (seconds)
    n_trials = cfg.n_trials;             % Number of trials
    t_trial = cfg.t_trial;               % Time vector for one trial
    f_osc = cfg.f_osc;                   % Oscillatory frequency (Hz)
    %event_time = cfg.event_time;         % Event times for each trial
    coupling_window = cfg.coupling_window; % Time window for coupling
    figure_flag = cfg.figure_flag;       % Flag to show figures
    spc_factor = cfg.spc_factor;         % Factor for spike-phase coupling
    phase_target = cfg.phase_target;     % Target phase for coupling
    fr_baseline = cfg.fr_baseline;       % Baseline firing rate
    gaussian_amplitude = cfg.gaussian_amplitude; % Gaussian modulation amplitude
    gaussian_sigma = cfg.gaussian_sigma; % Gaussian modulation standard deviation
    event_loc = cfg.event_loc;
    event_jitter = cfg.event_jitter;
    event2_loc = cfg.event2_loc;
    event2_jitter = cfg.event2_jitter;    
    % events
    cfg.event_time = event_loc + event_jitter * randn(1, cfg.n_trials); % Time of the event in each trial (seconds) + jitter
    cfg.event_time2 = event2_loc + event2_jitter * randn(1, cfg.n_trials); % Time of the other event in each trial (seconds)

    
    ES = (1 + 2*sqrt(spc_factor))./(1 - 2*sqrt(spc_factor));
    %ES = calc_effectsizePPC(spc_factor);
    cfg.ES = ES;
    % Preallocate arrays for results
    osc_signal_all_trials = zeros(n_trials, length(t_trial));
    spike_train_all_trials = zeros(n_trials, length(t_trial));
    time_windows = [0:0.05:(trial_duration-0.2); 0.2:0.05:trial_duration]'; 


    % Initialize random phase starts
    phase_start = -pi + 2 * pi * rand(1, n_trials);

    % Loop over trials
    for trial = 1:n_trials
        % Generate oscillatory signal (sinusoid)
        osc_signal = sin(2 * pi * f_osc * t_trial + phase_start(trial));
        osc_phase = angle(hilbert(osc_signal)); % Compute the phase

        % Baseline spike probability
        baseline_spike_prob = fr_baseline / fs; 

        % Gaussian modulation for transient effect
        gauss_window = (gaussian_amplitude - 1) * exp(-((t_trial - cfg.event_time(trial)).^2) / (2 * gaussian_sigma^2));
        
        % Create spike probability with modulation
        mod_spike_prob = baseline_spike_prob * (1 + gauss_window); 

        % Apply coupling only within the specified window
        %coupling_gaussian = exp(-((t_trial - mean(coupling_window)).^2) / (2 * (diff(coupling_window)/2.355)^2));
        coupling_idx = (t_trial >= coupling_window(1)) & (t_trial <= coupling_window(2));
        phase_difference = abs(osc_phase - phase_target);

        % Wrap phase difference to [0, pi]
        phase_difference(phase_difference > pi) = 2 * pi - phase_difference(phase_difference > pi);
        
        % Adjust spike probability based on phase difference
        phase_weighting = cos(phase_difference);
        
        % Redistribute spike probability within the coupling window
        spike_prob = mod_spike_prob;
        %spike_prob(coupling_idx) = spike_prob(coupling_idx) .* (1 + (ES - 1) * phase_weighting(coupling_idx));
        spike_prob(coupling_idx) = spike_prob(coupling_idx).*(1 + (ES - 1)/(ES + 1).*phase_weighting(coupling_idx));

        % Ensure non-negative spike probability
        %spike_prob = max(spike_prob, 0);

        % adjust mean
        %spike_prob = spike_prob - mean(spike_prob) + mean(mod_spike_prob);
        
        % Generate spike train based on probability
        spike_train = rand(size(t_trial)) < spike_prob;


        % Store results for this trial
        osc_signal_all_trials(trial, :) = osc_signal;
        spike_train_all_trials(trial, :) = spike_train;





    end

    % if align_toevent ~= 0
    %     t_trial_align = t_trial_align(logical(which_idx_trial(1,:)));
    %     osc_signal_all_trials = osc_signal_all_trials(:,logical(which_idx_trial));
    %     spike_train_all_trials = spike_train_all_trials(:,logical(which_idx_trial));
    % end

  

    % Calculate PLV and PPC across trials
    [PLV_across_trials, PPC_across_trials, nspikes] = calculate_PLV_traditional(spike_train_all_trials, osc_signal_all_trials, time_windows, t_trial, n_trials);
    time_bin = mean(time_windows, 2); % Calculate time bin centers
    SPC_method = calculate_PLV_PPC_win_simplified(spike_train_all_trials, osc_signal_all_trials, cfg);

    if cfg.nperms > 0
        SPC_method.Sign = SPC_method.PPC > prctile(squeeze(SPC_method.PPC_perm),95,2)';
    end
    % Calculate firing rate (FR)
    FR = nspikes / (0.2 * n_trials); 

    % Plot results if figure_flag is set
    if figure_flag
        plot_simulation(time_bin, FR, PPC_across_trials, SPC_method, cfg);
    end
end


