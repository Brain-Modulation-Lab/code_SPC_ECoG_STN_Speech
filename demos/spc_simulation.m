% Fixed random seed
rng(43)

% Define paths for script, data, and output
script_folder = pwd;  % Current folder of the script
[parentFolder,~,~] = fileparts(script_folder);


% Define paths for source scripts and data files
PATH_SPC = fullfile(parentFolder,'src');  % Path to source folder
addpath(genpath(PATH_SPC))

% Parameters
cfg = struct();
cfg.fs = 1000; % Sampling frequency (Hz)
cfg.trial_duration = 3; % Trial duration (seconds)
cfg.n_trials = 80; % Number of trials
cfg.t_trial = 0:1/cfg.fs:cfg.trial_duration - 1/cfg.fs; % Time vector for one trial
cfg.f_osc = 8; % Oscillatory frequency (Hz) (theta rhythm)
cfg.coupling_window = [1.3 1.7]; % Time window for spike-phase coupling (seconds)
cfg.spc_factor = 0.01; % Parameter to adjust occurrences based on phase
cfg.phase_target = pi/4;
cfg.fr_baseline = 25; % Baseline firing rate (Hz)
cfg.gaussian_amplitude = 1.5; % Peak firing rate (spikes per second) reduction by 50%
cfg.gaussian_sigma = 0.2; % Standard deviation for Gaussian (in seconds)
cfg.nperms = 0; %
cfg.window_bins = 0.025; % target size window (in ms)
cfg.figure_flag = false;
cfg.event_loc = 1.6;
cfg.event_jitter =  0;
cfg.event2_loc = 1.2;
cfg.event2_jitter =  0;

%
% Esperiment 0: Single run (SPC recovery)
cfg.figure_flag = true;
cfg.spc_factor = 0;
run_simulation(cfg)
cfg.spc_factor = 0.015;
run_simulation(cfg)
cfg.figure_flag = false;

%% Experiment 1: Effect of firing rate modulation
disp("Running experiment 1: Effect of firing rate modulation (with coupling):")
spc_factor = linspace(0,0.05,10);

% Run simulations with different firing rate modulations
cfg.gaussian_amplitude = 1; SPC1 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 0.5; SPC2 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 1.5; SPC3 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
%cfg.gaussian_amplitude = 0.5; SPC4 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);

% Plot results
SPCs = {SPC1, SPC2, SPC3};  % Collect all SPC results
plot_title = {"no mod","x0.5","x1.5"};
plot_experiment(SPCs, spc_factor,plot_title);

%%
figure("Position",[200 200 600 300])
tiledlayout(1,2)
nexttile
PPC = cellfun(@(x) arrayfun(@(y) ...
    mean(y.PPC_time(y.time_bin <= cfg.coupling_window(2) & ...
    y.time_bin >= cfg.coupling_window(1))), x), SPCs,"UniformOutput",false);
scatter([SPCs{1}.value],cell2mat(PPC')',75,'filled')
xlim([0 0.055])
ylim([0 0.055])
hold on
plot([0 0.0655],[0 0.055],'k--')
xlabel("True SPC")
ylabel("PPC")
legend(plot_title,'location','best')

PPCwin = nan(numel([SPCs{1}.value]),numel(SPCs));
for k = 1 : numel(SPCs)
    for ii = 1 : numel([SPCs{1}.value])
        timeVec = median(SPCs{k}(ii).WIN_METHOD.winCenters - SPCs{k}(ii).WIN_METHOD.winCenters(:, 1)) + median(SPCs{k}(ii).cfg.event_time2 - 0.8);
        PPC_vec = SPCs{k}(ii).WIN_METHOD.PPC;
        PPC_vec(PPC_vec<0) = 0;
        PPCwin(ii,k) = mean(PPC_vec(timeVec <= cfg.coupling_window(2) & timeVec >= cfg.coupling_window(1)));
    end
end


nexttile

scatter([SPCs{1}.value],PPCwin,75,'filled')
xlim([0 0.055])
ylim([0 0.055])
hold on
plot([0 0.0655],[0 0.055],'k--')
xlabel("True SPC")
ylabel("PPC win-method")

%plot_sweep_parameter(SPC)
%% Experiment 2: Effect of firing rate modulation w/ coupling elsewhere
disp("Running experiment 1: Effect of firing rate modulation (with coupling elsewhere):")
spc_factor = linspace(0,0.03,10);
ES = calc_effectsizePPC(spc_factor);

% Run simulations with different firing rate modulations
cfg.coupling_window = [1.1 1.4];
cfg.gaussian_amplitude = 1; SPC1 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 0.75; SPC2 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 1.25; SPC3 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 0.5; SPC4 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);

% Plot results
SPCs = {SPC1, SPC2, SPC3, SPC4};  % Collect all SPC results
cfg.plot_title = {"no mod","x0.75","x1.25","x0.5"};
plot_experiment(cfg, SPCs, spc_factor, ES);
%% Experiment 3: Effect of firing rate modulation w/ permutation
% example
spc_factor = 0.03;
ES = calc_effectsizePPC(spc_factor);

cfg.nperms = 1000;
cfg.gaussian_amplitude = 1; SPC1 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 0.75; SPC2 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 1.25; SPC3 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
cfg.gaussian_amplitude = 0.5; SPC4 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);

% Plot results
SPCs = {SPC1, SPC2, SPC3, SPC4};  % Collect all SPC results
cfg.plot_title = {"no mod","x0.75","x1.25","x0.5"};
plot_experiment(cfg, SPCs, spc_factor, ES);

%% Experiment 4: Effect of firing rate modulation w/ permutation
cfg.nperms = 500;
cfg.gaussian_amplitude = 1;
cfg.spc_factor = 0.015;

reps = 1;
SPC_sim = struct();

for rep = 1 : reps
    fprintf("Running repetition %d - %d \n",rep,reps)
    rng(41 + rep)
    % effect: ES
    spc_factor = linspace(0,0.04,10);
    %ES = (1 + 2 * sqrt(spc_factor)) ./ (1 - 2 * sqrt(spc_factor));
    SPC_ES2 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);
    %cfg.gaussian_amplitude = 1.25; SPC_ES3 = run_sweep_parameter(cfg, 'spc_factor', spc_factor);

    % effect: cluster duration
    coupling_window = 1.6 + [-[0.025:0.050 : 0.6]' [0.025:0.050 : 0.6]'];
    SPC_DUR2 = run_sweep_parameter(cfg, 'coupling_window', coupling_window);
    %cfg.gaussian_amplitude = 1.25; SPC_DUR3 = run_sweep_parameter(cfg, 'coupling_window', coupling_window);

    % effect: # bins
    window_bins = 0.01:0.02:0.1;
    SPC_BINS2 = run_sweep_parameter(cfg, 'window_bins', window_bins);
    %cfg.gaussian_amplitude = 1.25; SPC_BINS3 = run_sweep_parameter(cfg, 'window_bins', window_bins);

    % effect: window-width by changing trials
    n_trials = 60:10:90;
    SPC_WINWIDTH2 = run_sweep_parameter(cfg, 'n_trials', n_trials);
    %cfg.gaussian_amplitude = 1.1; SPC_WINWIDTH3 = run_sweep_parameter(cfg, 'n_trials', n_trials);

    % effect: Frequency oscillation
    f_osc = logspace(log10(4), log10(150), 10);
    %ES = (1 + 2 * sqrt(spc_factor)) ./ (1 - 2 * sqrt(spc_factor));
    SPC_FOSC2 = run_sweep_parameter(cfg, 'f_osc', f_osc);
    %vcfg.gaussian_amplitude = 1.25; SPC_FOSC3 = run_sweep_parameter(cfg, 'f_osc', f_osc);

    % effect: Average Firing rate
    fr_baseline = 5:5:60;
    SPC_FR2 = run_sweep_parameter(cfg, 'fr_baseline', fr_baseline);
    %cfg.gaussian_amplitude = 1.25; SPC_FR3 = run_sweep_parameter(cfg, 'fr_baseline', fr_baseline);


    % effect: jitter events2
    event_jitter = 0:0.025:0.15;
    SPC_EVTJIT2 = run_sweep_parameter(cfg, 'event_jitter', event_jitter);
    %cfg.gaussian_amplitude = 1.25; SPC_EVTJIT3 = run_sweep_parameter(cfg, 'event_jitter', event_jitter);


    % effect: firing modulation
    gaussian_amplitude = 0.25:0.25:4;
    SPC_FRMOD2 = run_sweep_parameter(cfg, 'gaussian_amplitude',gaussian_amplitude);

    % effect: fs
    fs = 200:200:1000;
    SPC_FS2 = run_sweep_parameter(cfg, 'fs', fs);
    %cfg.gaussian_amplitude = 1.25; SPC_FS3 = run_sweep_parameter(cfg, 'fs', fs);


    % store reps
    SPC_sim(rep).id = rep;
    SPC_sim(rep).spc_factor.store = SPC_ES2;
    SPC_sim(rep).coupling_window.store = SPC_DUR2;
    SPC_sim(rep).window_bins.store = SPC_BINS2;
    SPC_sim(rep).n_trials.store = SPC_WINWIDTH2;
    SPC_sim(rep).f_osc.store = SPC_FOSC2;
    SPC_sim(rep).fr_baseline.store = SPC_FR2;
    SPC_sim(rep).event_jitter.store = SPC_EVTJIT2;
    SPC_sim(rep).gaussian_amplitude.store = SPC_FRMOD2;
    SPC_sim(rep).fs.store = SPC_FS2;
end

%% plot results experiment 4
params = ["spc_factor","coupling_window","window_bins","n_trials","f_osc","fr_baseline","event_jitter","gaussian_amplitude","fs"];
params_label = ["true PPC","PPC window [s]","# bins","# trials","Freq. osc. [Hz]","Firing rate [Hz]","Jitter [s]","Firing rate mod. ","Sampling rate [Hz]"];


figure
for p= 1 : numel(params)
    nexttile
    hold on
    param = params(p);
    for rep = 1 : reps
        SPC_sim_ = SPC_sim(rep).(param);
        if strcmpi(param,"coupling_window")
            value = arrayfun(@(x) diff(x.value),SPC_sim.(param).store);
            range_par = [0.95*min(value) 1.05*max(value)];

        else
            value = [SPC_sim.(param).store.value];
            range_par = [0.95*min(value) 1.05*max(value)];
        end
        if ~strcmpi(param,"window_bins") && ~strcmpi(param,"n_trials")
            PPC_mat = [];
            for jj = 1 : numel(SPC_sim_.store)
                SPCwin = SPC_sim_.store(jj).WIN_METHOD;
                PPC_mat = [PPC_mat SPCwin.PPC'];
                timeVec = median(SPCwin.winCenters - SPCwin.winCenters(:, 1)) + mean(SPCwin.cfg.event_time2 - 0.8);

            end
            pcolor(value,timeVec,PPC_mat/6);
            colormap(linspecer)
            cb = colorbar;
            ylabel(cb,"PPC [raw]")
        end
        for jj = 1 : numel(SPC_sim_.store)
            SPCwin = SPC_sim_.store(jj).WIN_METHOD;
            timeVec = median(SPCwin.winCenters - SPCwin.winCenters(:, 1)) + mean(SPCwin.cfg.event_time2 - 0.8);
            SPCwin.Sign = SPCwin.PPC > prctile(squeeze(SPCwin.PPC_perm),95,2)';

            scatter(value(jj)*ones(1,sum(SPCwin.Sign)),  timeVec(SPCwin.Sign),'SizeData',75,'Marker','x','MarkerEdgeColor','k','LineWidth',2)
            store = comp_clusterPLV(SPC_sim_.store(jj).WIN_METHOD.PPC, SPC_sim_.store(jj).WIN_METHOD.PPC_perm, 'zstat');
            scatter(value(jj)*ones(1,sum(store.LabelMatrix>0)),  timeVec(store.LabelMatrix>0),'SizeData',75,'Marker','x','MarkerEdgeColor','r','LineWidth',2)
        end
    end
    if strcmpi(param,"coupling_window")
        patch([value(1) value(end) value(end)], 1.6 + [0  -max(value)/2   max(value)/2], 'r', 'FaceAlpha', 0.3);
    else
        patch([range_par fliplr(range_par)], [cfg.coupling_window(1) cfg.coupling_window(1) cfg.coupling_window(2) cfg.coupling_window(2) ], 'r', 'FaceAlpha', 0.3);
    end
    %scatter(value,  squeeze(Toff(:,1,:))',75,'filled')
    xlabel(params_label{p},"Interpreter","none")
    ylabel("SPC time [s]")
    ylim([1 2.3])
    xticks(round(value(1:2:end),2))
end

