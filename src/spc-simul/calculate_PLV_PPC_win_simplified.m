function store = calculate_PLV_PPC_win_simplified(spike_train_all_trials, osc_signal_all_trials, cfg)
    % Function to calculate Phase Locking Value (PLV) and Phase Phase Coherence (PPC) for given spike trains and oscillatory signals.

    ntrials = cfg.n_trials; % Number of trials
    spike_train_all_trials = spike_train_all_trials'; % Transpose spike train matrix
    spike_train_all_trials = spike_train_all_trials(:)'; % Reshape to 1D array

    % Extract epochs of spike data aligned to main event of interest
    [spike_atSpeech_trials] = getTrialsMatrix(cfg.event_time' + cfg.trial_duration * [0:(cfg.n_trials-1)]', ...
                                                spike_train_all_trials, 1/cfg.fs, 2, 1); 
    spikeTrialSum = sum(spike_atSpeech_trials, 2); % Sum spikes across trials
    nSpikes = round(sum(spikeTrialSum) / 2 * 0.2); % Calculate number of spikes for analysis

    %% Compute the time points of interest for each trial
    allEvts = []; % Initialize event matrix
    % Define event time points with offsets for windows
    allEvts = [cfg.event_time2' - 0.8, cfg.event_time2', cfg.event_time', cfg.event_time' + 0.6] + ...
              cfg.trial_duration * [0:(cfg.n_trials-1)]'; 
    allEvts = allEvts'; % Transpose event matrix

    % Decide number of bins for analysis
    num_bins = floor(mean(cfg.event_time + 0.6 - (cfg.event_time2 - 0.8)) / cfg.window_bins / 4) * ones(3, 1);
    winCenters = []; % Initialize window centers

    % Loop through event time points to generate window centers
    for e = 1:(size(allEvts, 1) - 1)
        winCenters_tmp = []; % Temporary storage for window centers
        for trl = 1:size(allEvts, 2)
            % Create equidistant points between neighboring events
            winCenters_tmp(trl, :) = linspace(allEvts(e, trl), allEvts(e + 1, trl), num_bins(e));
        end
        if e == 1
            winCenters = [winCenters, winCenters_tmp]; % Append first set of centers
        else
            % Remove first point of current centers to avoid duplication with last point of previous centers
            winCenters = [winCenters, winCenters_tmp(:, 2:end)];
        end
    end

    %% Extract the LFP phase for the frequency of interest
    phase = angle(hilbert(osc_signal_all_trials))'; % Compute phase of the oscillatory signal using Hilbert transform
    phase = phase(:)'; % Reshape phase to 1D array

    % Handle power trimming (commented out, but could be relevant)
    spike_Vect_foi = spike_train_all_trials; % Store original spike train vector

    %% Calculate PLV in the windows and estimate the bin length for each window
    idcs = round(winCenters(2:end, :) * cfg.fs); % Convert window centers to sample indices
    store = struct(); % Initialize output structure

    % Loop through each window center
    for w = 1:size(winCenters, 2)
        % Function to select spikes based on bin length and compute spikes
        [curr_spks, store] = calc_binLen_selectSpikes_simplified(idcs, w, spike_Vect_foi, nSpikes, store, nan, nan);

        currPhases = phase(curr_spks == 1); % Extract phases corresponding to current spikes
        [PLV_original, ~] = comp_PLV(currPhases); % Compute original PLV
        store.PLV(w) = PLV_original; % Store PLV result
        store.PPC(w) = correct_PLV(PLV_original, sum(~isnan(currPhases))); % Store corrected PLV (PPC)
        if w == 31
            figure
            histogram(currPhases)
        end
        % Perform permutation procedure for significance testing
        [phases_atSpeech_trls, phase_atSpeech_vect] = getTrialsMatrix(winCenters(2:end, w), phase, 1/cfg.fs, 1, 0.5);
        
        for perm = 1:cfg.nperms % Permutation test loop
            rng(perm); % Set random seed for reproducibility
            randIdcs = randperm(size(phases_atSpeech_trls, 2)); % Generate random indices for shuffling
            shuffledPhaseTrials = phases_atSpeech_trls(:, randIdcs); % Shuffle phase trials
            tmp = phase_atSpeech_vect; % Store original phase vector
            tmp(~isnan(tmp)) = shuffledPhaseTrials(:); % Replace with shuffled trials
            shuffledPhase = tmp; % Create shuffled phase vector
            [PLV_perm, ~] = comp_PLV(shuffledPhase(curr_spks == 1)); % Compute PLV for shuffled phases
            store.PLV_perm(w, perm) = PLV_perm; % Store permutation PLV
            store.PPC_perm(w, perm) = correct_PLV(PLV_perm, sum(~isnan(shuffledPhase(curr_spks == 1)))); % Store corrected PPC for shuffled
        end
    end

    % Save global structure for output
    store.cfg = cfg; % Store configuration
    store.n_trials = ntrials; % Store number of trials
    store.winCenters = winCenters; % Store calculated window centers
    store.centerEvts = cumsum([1; num_bins - 1]); % Cumulative sum for x-Ticks in plots

end
