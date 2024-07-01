%% This code reproduces the Figure 2 of the manuscript. This figure describes
% the properties of the identified events of spike-phase coupling.

clc
clear all
close all

% Define paths for script, data, and output
script_folder = pwd;  % Current folder of the script
[parentFolder,~,~] = fileparts(script_folder);


% Define paths for source scripts and data files
PATH_SPC = fullfile(parentFolder,'src');  % Path to source folder
PATH_DATA = fullfile(parentFolder,'data');  % Path to source folder
PATH_DB = fullfile(PATH_DATA,'DB_main_analysis.mat');  % Path to source folder

% Add source path and its subdirectories to the search path
addpath(genpath(PATH_SPC))

% initialize libraries
bml_defaults;
ft_defaults;

% make plot prettier
prettify_plots;

% define frequency bands
bands = bml_get_canonical_bands([0,150]);
bands.fmid = sqrt(bands.fstarts .* bands.fends);
newcolor = {'#d7191c', '#fdae61','#fecc5c','#abdda4','#2b83ba'};
bands.color(2:6) = newcolor;


% Define your participant list
SUBJECTS = {'DBS3001','DBS3002','DBS3003','DBS3004','DBS3005','DBS3008', 'DBS3010', 'DBS3011', ...
            'DBS3012', 'DBS3014', 'DBS3015', 'DBS3017', 'DBS3018', 'DBS3019', 'DBS3020', 'DBS3022', ...
            'DBS3023', 'DBS3024', 'DBS3026', 'DBS3027', 'DBS3028', 'DBS3029', 'DBS3030', 'DBS3031', ...
            'DBS3032'};


% Step 1. gather all SPC (pairs & clusters) information
load(PATH_DB,'DB')
cfg = set_configs('default');


% Step 2. Get times of events
EvtTimes = get_EvtTimes(DB.Pairs);

% Step 3. Extract SPC maps ('This takes a while')
Evt_lock = 'Speech'; % or Cue
which_pairs = 'OnlySign'; % or 'OnlySign', 'Top'
which_var = 'PPCz';
Tres = 0.05; % resolution of map
T = [-2.5 2]; % rnage of time interval
nPerms = 500; % number of permutation

cfg.stat_test = 'baseline-test'; % can be 'paired-test-conds' (use this if you want to comapre two conditions), 'baseline-test'
cfg.Filestats = fullfile(PATH_DATA,'permutation_avgmaps',[which_var,'Stat_pairs-',which_pairs,'_evt-',Evt_lock,'_nperm-',num2str(nPerms),'.mat']); % use the pre=saved permutaiton because they take long time!
%cfg.Baseline = [-Inf EvtTimes.Speech(2)];

SPCmap = get_SPCmap(DB.Pairs,which_var, T, ...
    Tres,cfg,Evt_lock,which_pairs,nPerms);

% Panel Figure 2A: Plot average SPC maps 
cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTimes = EvtTimes.Speech;
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = {[0 0.7],[-10 10]};
cfg_fh.CMap = linspecer;
cfg_fh.NClustplot = 4;
fhA = plotter_SPCmap(SPCmap, cfg_fh);

% Panel Figure 2B-C: Plot all single pairs in a subject [DBS3008]
% take the pairs
cfg_fh = struct();
cfg_fh.freqsOfInt = cfg.plv.freqsOfInt;
cfg_fh.Subjects = SUBJECTS(strcmpi(SUBJECTS,'DBS3008')); % 
SPCpair = get_SPCpair(DB, cfg_fh); % get only significant pairs
% plot them
cfg_fh = struct();
cfg_fh.Position = [300 300 700 500];
cfg_fh.EvtTypes = cfg.plot.EvtTypes;
cfg_fh.CLim = [-2.5 2.5]; % hard to set a priori
cfg_fh.CMap = linspecer;
cfg_fh.Visible = 'off';
cfg_fh.smooth = true; % gaussian smoothing to improve visualization
fhB = plotter_SPCpair(SPCpair, cfg_fh); % cell of figures handles

% Panels Figure 2D-F: Plot t-SPC duration and frequency

clust_res = .005;
Tc = -2.5 : clust_res : 2;

% Panel D
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = Tc;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'Duration_1DHist','Frequency_1DHist'};
fhD = plot_ClusterProperties(cfg_cl, DB.Clusters);

% Panel E
cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = Tc;
cfg_cl.VarFields = {'ClusterCentroid','ClusterCentroidBand','ClusterOnOff'};
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.unit_type = 'all';
fhE = plot_ClusterTime(cfg_cl,DB.Clusters);

% Panel F (this is slow, it can take time (~15 min on personal laptop, ~3-5 min on our server)!)

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = Tc;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'TaskModulation_1DAggregationLine'};
fhF = plot_ClusterProperties(cfg_cl, DB.Clusters);

% Panel F (alternative with stacked areas, way faster!)

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = Tc;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.Properties = {'TaskModulation_1DArea'};
fhF_v2 = plot_ClusterProperties(cfg_cl, DB.Clusters);

% Panel Figure G: Plot frequency specificity
cfg_cl = [];
cfg_cl.bands = bands;
cfg_cl.freqsOfInt = cfg.plv.freqsOfInt;
cfg_cl.nPerms = nPerms;
fhG = plot_FrequencySpecificity(cfg_cl, DB);

% Supplementary Figure S7: Time delay

cfg_cl = [];
cfg_cl.cfg = cfg;
cfg_cl.Task_labels = {'ITI','Cue','Pre-Speech','Speech','Post-Speech'};
cfg_cl.T = Tc;
cfg_cl.EvtTimes = EvtTimes.Speech;
cfg_cl.bands = bands;
cfg_cl.freq_range = [4 30];
cfg_cl.Properties = {'TaskModulation_LagDelay','TaskModulation_Phasepolar','TaskModulation_Lagdelay_overtime'};
fhS7 = plot_ClusterProperties(cfg_cl, DB.Clusters);


