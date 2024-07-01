function store_all = calc_spike_PLV_all(E, S, cfg_analysis)

WIN_GET_SPIKENR = cfg_analysis.plv.WIN_GET_SPIKENR; % [sec], if you do not have lots of trials and thus not so many spikes you will want to make this bigger
tWidth_avgSpike = cfg_analysis.plv.tWidth_avgSpike;  % [sec] window size and..
tOffset_avgSpike = cfg_analysis.plv.tOffset_avgSpike;  % [sec] ...offset used to calculate the average spike number

NUM_PERMS = cfg_analysis.plv.NUM_PERMS;       % number of permutations, set at least to 200
PERM_WIN_LEN_SEC = cfg_analysis.plv.PERM_WIN_LEN_SEC; % this is used to permute the phase-signal, the exact number does not matter,
% it only needs to be longer than the largest bin that you'll use to compute the PLV with


LENGTH_WINDOW = cfg_analysis.plv.LENGTH_WINDOW;  % this is to set NBINS_WINDOW of bins in a window LENGTH_WINDOW
NBINS_WINDOW = cfg_analysis.plv.NBINS_WINDOW;
THR_NAN = cfg_analysis.plv.THR_NAN;
MIN_TRIALS = cfg_analysis.plv.MIN_TRIALS;
MAX_NAN = cfg_analysis.plv.MAX_NAN;

freqsOfInt = cfg_analysis.plv.freqsOfInt;   % Define your frequencies of interest

% find channel labels
E_channels = E.label;
S_channels = S.label;
% find all possible pairs
[tempA, tempB] = ndgrid(1 : numel(E_channels) , 1: numel(S_channels));
pairs_idx = [tempA(:) tempB(:)];
npairs = size(pairs_idx,1);

% LFP_channel   = 'LFP1';  % rename those according to your data
% spike_channel = 'unit1';

% loop through pairs
store_all = struct();

count_pair = 1;
for pair_i = 1:npairs
    
    clear store
    store = struct();
    
    % take channels and store them
    E_channel = E_channels{pairs_idx(pair_i,1)};
    S_channel = S_channels{pairs_idx(pair_i,2)};
    
    
    % select neuron
    idx_neuron = find(strcmpi(S.label,S_channel));
    spikeTimes = S.timestamp{idx_neuron};
    % need coloumn vector for spike times
    if size(spikeTimes,1) == 1
        spikeTimes = spikeTimes';
    end
    
    % spike times in seconds
    if ~any(strcmp('session_id_coding',E.epochs.Properties.VariableNames))
       E.epochs.session_id_coding = E.epochs.session_id;
    end
     
    try
        sess_neuron = unique(S.epochs{idx_neuron}.session_id_coding);
    catch
        sess_neuron = unique(S.epochs{idx_neuron}.session_id);
        S.epochs{idx_neuron}.session_id_coding = S.epochs{idx_neuron}.session_id;        
    end
    % just a check
    assert(numel(sess_neuron) == 1, "Single-Unit is present during more ECoG sessions - need to handtune the code!")
    
    
    
    % select E channel  during proper session
    cfg_ft = [];
    cfg_ft.channel = {E_channel};
    E_data = ft_selectdata(cfg_ft,E);
    
    % correct if locked
    if cfg_analysis.locked
        E_epochs = E.epochs(E.epochs.session_id_coding == sess_neuron,:);
        x_ = [];
        for ee = 1 : height(E_epochs)
            start_ = E_epochs.starts(ee);
            stop_ = E_epochs.ends(ee);
            idx_ = E_data.time{sess_neuron} >= start_ & E.time{sess_neuron} <= stop_;
            x_ = [x_ ;E_data.trial{sess_neuron}(idx_)];
        end
        evoked = mean(x_,'omitnan');
        for ee = 1 : height(E_epochs)
            start_ = E_epochs.starts(ee);
            stop_ = E_epochs.ends(ee);
            idx_ = E_data.time{sess_neuron} >= start_ & E.time{sess_neuron} <= stop_;
            E_data.trial{sess_neuron}(idx_) = E_data.trial{sess_neuron}(idx_) - evoked;
        end
    end
    
    
    % select current data session
    E_data = E_data.trial{sess_neuron};
    
    
    
    nan_perc = 1E2*sum(isnan(E_data))/numel(E_data);
    % select sampling rate
    SR = E.fsample;
    
    % get neurons events (need to correct for non-nan ecog epochs, it could be brutal the correction losing so many trials...)
    epochs = S.epochs{idx_neuron};
    
    % find epochs with nan in E_data
    time_nan = E.time{sess_neuron}(isnan(E_data));
    try
        flag_nan = sum(time_nan >= epochs.starts_coding & time_nan <= epochs.ends_coding,2) > THR_NAN;
    catch
        flag_nan = sum(time_nan >= epochs.starts & time_nan <= epochs.ends,2) > THR_NAN;
    end
    epochs_valid = epochs(~flag_nan,:);
    
    % if spikes lands in
    ntrials = height(epochs_valid);
    
    
    
    
    if ntrials >= MIN_TRIALS && nan_perc <= MAX_NAN
        fprintf('Pair %d: %s-%s # trials: %i\n',pair_i, E_channel, S_channel, ntrials) %currFile, LFP_channel, spike_channel, freqOfInt(1));
        
        % adjust spike-times and events by first E timing
        offTime = E.time{sess_neuron}(1);
        spikeTimes = spikeTimes - offTime;
        
        Evt = table();
        Evt.cue_onset = epochs_valid.stim1_onset - offTime;
        Evt.cue_offset = epochs_valid.stim3_offset - offTime;
        Evt.speech_onset = epochs_valid.syl1_onset - offTime;
        Evt.speech_offset = epochs_valid.syl3_offset - offTime;
        
        %     currFile = fileNames{f};
        %     load([Paths.DataDir, currFile,'.mat'])
        %     Evt = data.Evt;
        %     LFP = data.LFP;
        %     SR  = data.LFP_SR;
        %     spikeTimes = data.spikeTimes;  % spike times in seconds
        
        %% Load spikes
        spike_Vect = zeros(1, length(E_data));
        idcs_tmp = round(spikeTimes*SR);
        idcs_tmp(idcs_tmp == 0) = []; % remove the index of 0 if there was one
        spike_Vect(idcs_tmp) = 1;
        
        %% Get epochs around your main event of interest to calculate the number of spikes you want to include per window
        % based on the average within a window of a duration defined in WIN_GET_SPIKENR
        [spike_atSpeech_trials]  =  getTrialsMatrix(Evt.cue_onset, spike_Vect, 1/SR, tWidth_avgSpike, tOffset_avgSpike);  % get epochs, aligned to your main event of interest
        spikeTrialSum = sum(spike_atSpeech_trials,2);
        nSpikes       = round(sum(spikeTrialSum) / tWidth_avgSpike * WIN_GET_SPIKENR);
        
        
        %% Compute the time points of interest for each trial
        allEvts = [];
        trls_speech = ~isnan(Evt.cue_onset + Evt.cue_offset + Evt.speech_onset + Evt.speech_offset);
        %trls_mvmt(end) = false; % remove the last event because the data was too short in some cases, @Matteo, you can delete this line
        
        allEvts(1,:) = Evt.cue_onset(trls_speech) - 0.75;
        allEvts(2,:) = Evt.cue_onset(trls_speech);
        allEvts(3,:) = Evt.cue_offset(trls_speech);
        allEvts(4,:) = Evt.speech_onset(trls_speech);
        allEvts(5,:) = Evt.speech_offset(trls_speech);
        allEvts(6,:) = Evt.speech_offset(trls_speech) + 0.75;
        
        num_bins = round(diff(mean(allEvts,2)) / LENGTH_WINDOW * NBINS_WINDOW); % decide how many bins you want to have, this here is set to create 10 bins per 0.5 seconds
        
        winCenters = [];
        for e = 1:(size(allEvts,1)-1)
            winCenters_tmp = [];
            for trl = 1:size(allEvts,2)
                % Get the appropriate centers for each trial by
                % making equidistant points between neighbouring
                % events
                winCenters_tmp(trl,:) = linspace(allEvts(e,trl), allEvts(e+1,trl), num_bins(e));
            end
            if e == 1
                winCenters = [winCenters, winCenters_tmp]; % simply append the first time
            else % and for any subsequent events, remove the first winCenter point, as this is
                % identical to the last winCenter point of the previous event
                winCenters = [winCenters, winCenters_tmp(:, 2:end)];
            end
        end


        if cfg_analysis.flip.flip_check
            flip_flag = flip_polarity(E_data, SR, cfg_analysis.flip);

        end

        if cfg_analysis.plv.gamma_env
            [bFOI, aFOI] = butter(cfg_analysis.flip.ordFfoi, cfg_analysis.plv.gamma_env_BPF/(SR/2));
            tmp= E_data;
            tmp(isnan(E_data)) = 0;
            GammaLFP_filt = filtfilt(bFOI, aFOI,tmp) ;
            GammaLFP_filt = abs(hilbert(GammaLFP_filt));
        end


        for frq = 1:size(freqsOfInt,1)
            freqOfInt = freqsOfInt(frq,:);
            
            %% Extract the LFP phase for the frequency of interest
            %         LFP_filt1 = ft_preproc_highpassfilter(LFP, SR, freqOfInt(1), 4, 'but','twopass'); % 3rd order butterworth high-pass filter
            %         LFP_filt1 = ft_preproc_lowpassfilter(LFP_filt1, SR, freqOfInt(2), 4, 'but','twopass'); % 3rd order butterworth low-pass filter
            %         d = designfilt('bandpassiir', 'FilterOrder',4, ...
            %             'HalfPowerFrequency1',freqOfInt(1),'HalfPowerFrequency2',freqOfInt(2), ...
            %             'SampleRate',SR);
            %         LFP_filt = filter(d, LFP);  %% !!! change this into a two-pass filter (passed forwards and backwards, so that you'll get zero phase distortion)
            %         LFP_filt = filtfilt(d, LFP);  % watch out, this gives weird
            %         results if you have the fieldtrip toolbox or EEGLAB toolbox added
            %         to your path, currently not working
            
            if ~cfg_analysis.plv.gamma_env
                [bFOI, aFOI] = butter(cfg_analysis.flip.ordFfoi, freqOfInt/(SR/2));
                tmp= E_data;
                tmp(isnan(E_data)) = 0;
                LFP_filt = filtfilt(bFOI, aFOI,tmp) ;
                phase = angle(hilbert(LFP_filt));
            else
                [bFOI, aFOI] = butter(cfg_analysis.flip.ordFfoi, freqOfInt/(SR/2));
                tmp = GammaLFP_filt;
                tmp(isnan(GammaLFP_filt)) = 0;
                LFP_filt = filtfilt(bFOI, aFOI,tmp) ;
                phase = angle(hilbert(LFP_filt));                
            end

            % handle power trim
            if  cfg_analysis.powertrim
                envtmp = abs(hilbert(LFP_filt));
                switch  cfg_analysis.powertrim_direction
                    case 'low'
                envtmp_low = envtmp <= prctile(envtmp(~isnan(E_data)),cfg_analysis.powertrim_prctile);
                    case 'high'
                  envtmp_low = envtmp >= prctile(envtmp(~isnan(E_data)),cfg_analysis.powertrim_prctile);              
                end
                    envtmp_low = setdiff(find(envtmp_low),find(isnan(E_data)));
                spike_Vect_foi = spike_Vect;
                spike_Vect_foi(envtmp_low) = 0;
            else
                spike_Vect_foi = spike_Vect;
            end

            if flip_flag
                phase = calc_shiftPhase(phase,pi);
            end
            %phase(isnan(E_data)) = nan;
            % check phase correct
            
            
            %% Calculate PLV in the windows and estimate the binLength (or window width) for each window
            idcs = round(winCenters*SR);
            fprintf('Pair %d: %s-%s @freq: %i\n',pair_i, E_channel, S_channel, freqOfInt(1)) %currFile, LFP_channel, spike_channel, freqOfInt(1));
            for w = 1:size(winCenters,2)
                % Important function that finds the appropriate
                % window-width to match the number of spikes across
                % windows that are used to calculate the PLV
                [curr_spks, store] = calc_binLen_selectSpikes(idcs, w, spike_Vect_foi, nSpikes, store, nan, nan);
                
                currPhases                = phase(curr_spks==1);
                [PLV_original, meanAngle] = comp_PLV(currPhases);
                store.PLV(frq, w)         = PLV_original;
                store.phase(frq, w)       = meanAngle;
                store.PPC(frq, w)          = correct_PLV(PLV_original, sum(~isnan(currPhases))); % correct plv
                store.ES(frq, w)           = calc_effectsizePPC(store.PPC(frq, w)); % correct effect size
                
                % Get the phase for the permutation procedure:
                %phase(isnan(phase)) = 0;
                [phases_atSpeech_trls, phase_atSpeech_vect]  = getTrialsMatrix(winCenters(:,w), phase, 1/SR,PERM_WIN_LEN_SEC, PERM_WIN_LEN_SEC/2);
                for perm = 1:NUM_PERMS
                    rng(perm);  % use rng here, otherwise the perm distribution
                    % would be re-shuffled for each window/frequency
                    randIdcs = randperm(size(phases_atSpeech_trls,2));  % create random indices to shuffle the order of trials
                    shuffledPhaseTrials = phases_atSpeech_trls(:,randIdcs);
                    tmp              = phase_atSpeech_vect;       % this is the full time series, but only includes
                    tmp(~isnan(tmp)) = shuffledPhaseTrials(:);  % refill the long vector with the shuffled trials
                    shuffledPhase    = tmp;
                    [PLV_perm, meanAngle_perm]     = comp_PLV(shuffledPhase(curr_spks==1));
                    store.PLV_perm(frq, w, perm)   = PLV_perm;
                    store.phase_perm(frq, w, perm) = meanAngle_perm;
                    store.PPC_perm(frq, w, perm)   = correct_PLV(PLV_perm, sum(~isnan(shuffledPhase(curr_spks==1))));
                end
            end
        end
        
        % Save the global structure
        store_all(count_pair).id = pair_i;
        store_all(count_pair).cfg = cfg_analysis;
        store_all(count_pair).E_channel = E_channel;
        store_all(count_pair).S_channel = S_channel;
        store_all(count_pair).n_trials = ntrials;
        store_all(count_pair).winCenters = winCenters;
        store_all(count_pair).centerEvts = cumsum([1; num_bins-1]); % this will later be used to set the x-Ticks and x-Tick labels in the time-frequency plots
        store_all(count_pair).freqsOfInt   = freqsOfInt;
        if cfg_analysis.flip.flip_check
            store_all(count_pair).flip_polarity = flip_flag;
        else
            store_all(count_pair).flip_polarity = nan;
        end
        store_all(count_pair).store = store;
        
        count_pair = count_pair + 1;
        
    else
        fprintf('Skipped Pair %d: %s-%s # trials: %i\n',pair_i, E_channel, S_channel, ntrials) %currFile, LFP_channel, spike_channel, freqOfInt(1));
        %         store_all(pair_i).id = pair_i;
        %         store_all(pair_i).cfg_analysis= cfg;
        %         store_all(pair_i).E_channel = E_channel;
        %         store_all(pair_i).S_channel = S_channel;
        %         store_all(pair_i).n_trials = ntrials;
        %         store_all(pair_i).winCenters = nan;
        %         store_all(pair_i).centerEvts = nan; % this will later be used to set the x-Ticks and x-Tick labels in the time-frequency plots
        %         store_all(pair_i).freqsOfInt   = nan;
        %         store_all(pair_i).flip_polarity = nan;
        %         store_all(pair_i).store = nan;
    end
    
    % #########################
    %     %% Save the file
    %     if ~exist([Paths.saveMatFiles, '/PLV'])
    %         mkdir([Paths.saveMatFiles, '/PLV'])
    %     end
    %     save([Paths.saveMatFiles, '\PLV\PLV_', currFile, '_', LFP_channel, '_', spike_channel, 'NUM_PERM=', num2str(NUM_PERMS), '_WIN=', sprintf('%.2f', WIN_GET_SPIKENR), '.mat'],  'store', '-v7.3')
end
